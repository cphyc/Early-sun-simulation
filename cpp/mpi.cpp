#include "mpi.h"
#include <iostream>
#include "RunnerModule.hpp"
#include "simul_constants.hpp"

inline int min(int x, int y) {
  return (x < y)? x : y;
}
/* Manager code create the data and send it one by one to 
   the workers. In the end, it fetches the result, print it to
   the output and ask to exit */
void manager_code(int numprocs);

/* Take no input. Wait for a message through MPI (a row of dOmegas), then process
   it, return the answer and wait for a new task. */
void worker_code();

/* Just print the FGM */
void print_FGM(double FGM[NOMEGA][NDOMEGA]);


/* Calculte the ranges */
void set_ranges() {
  // Set the ranges
  // logarithmic scale for k
  double kmin = 2*simul::pi/simul::lmax;
  double kmax = 2*simul::pi/simul::lmin;
  double fact = pow(kmax/kmin, 1./(NK-1));
  for (int n = 0; n < NK; n++) {
    double tmp = kmin * pow(fact, n);
    simul::k_range[n]    = tmp;
    simul::k_range[NK+n] = -tmp;
  }

  // Omega range and dlnOmegadlnr_range
  for (int n = 0; n < NOMEGA; n++)
    simul::Omega_range[n] = 31./NOMEGA*simul::OmegaSun*n;
  for (int n = 0; n < NDOMEGA; n++)
    simul::dlnOmegadlnr_range[n] = -2.5/NDOMEGA*n;
}
int main(int argc, char *argv[]) {
  int numprocs, rank;

  set_ranges();

  MPI_Init(&argc, &argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if (rank == 0) // Manager
    manager_code(numprocs);
  else
    worker_code();
  MPI_Finalize();
  return 0;
}

void manager_code(int numprocs){
  MPI_Status status;
  int n, nsent = 0;
  int om_rank, sender;
  double ret[NDOMEGA];
  
  // Creates the input array ...
  double domegas[NDOMEGA+1];

  for (int i = 0; i < NDOMEGA; i++)
    domegas[i] = simul::dlnOmegadlnr_range[i];

  // ... and the output array !
  double FGM[NOMEGA][NDOMEGA];

  /* Send a first row (at most numprocs, or NOMEGA if too many numprocs
     np : the tag 0 is used to send the exit message */
  for ( int i = 1; i <= min(numprocs, NOMEGA); i++) {
    // build the argument (don't forget, i starts a 1)
    domegas[NDOMEGA] = simul::Omega_range[i-1];
    // Send a column of domega at constant omega (last element : omega)
    MPI_Send(domegas, NDOMEGA + 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
    nsent++;
  }

  // Receive the results (starting with previous n !)
  for ( int i = 0; i < NOMEGA; n++) {
    /* Get back the answer omega rank via status.MPI_TAG and a sender number
       from any source. */
    MPI_Recv( &ret, NDOMEGA, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG,
	      MPI_COMM_WORLD, &status);
    /* Fetch back the omega rank and the sender.*/
    sender  = status.MPI_SOURCE;
    om_rank = status.MPI_TAG - 1;

    /* Store the answer in our array */
    for (int dom_rank = 0; dom_rank < NDOMEGA; dom_rank++) {
      FGM[om_rank][dom_rank] = ret[dom_rank];
    }

    /* Send another work (if any) ! */
    if (nsent < NOMEGA){
      domegas[NDOMEGA] = simul::Omega_range[nsent];
      MPI_Send(domegas, NDOMEGA + 1, MPI_DOUBLE, n, n, MPI_COMM_WORLD);
      nsent ++;
    }
    /* No other works, tell the sender to stop */
    else
      MPI_Send(MPI_BOTTOM, 0, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD);
  }
  print_FGM(FGM);
}


void worker_code(void) {
  MPI_Status status;
  int rank, Omega_rank;
  double Omega, dOmega;
  double FGM[NDOMEGA];
  double in[NDOMEGA+1];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  /* Should always be true ... */
  if (rank <= NOMEGA) {
    /* Get a new command from tag 0 with a null (0) tag */
    MPI_Recv(in, NOMEGA+1, MPI_DOUBLE, 0, 0, 
	     MPI_COMM_WORLD, &status);
    
    while (status.MPI_TAG > 0) {
      Omega_rank = status.MPI_TAG - 1;
      // Get the last element of the command as omega
      Omega = in[NDOMEGA];
      
      /* Call the get_FGM function that finds the fastest 
	 growing mode at a given Omega by looping over dOmega */
      for (int i; i < NDOMEGA; i++) {
	dOmega = in[i];
	FGM[i] = get_FGM(Omega, dOmega);
      }
      
      /* Send the answer (42) to 0 with the Omega_rank.*/
      MPI_Send(&FGM, NDOMEGA, MPI_DOUBLE, 0, Omega_rank + 1,
	       MPI_COMM_WORLD);

      /* Get a new job */
      MPI_Recv(&in, NDOMEGA+1, MPI_DOUBLE, 0, 0, 
	       MPI_COMM_WORLD, &status);
    }
  }
}

void print_FGM(double FGM[NOMEGA][NDOMEGA]){
  for (int i = 0; i < NOMEGA; i++) {
    for (int j = 0; j < NDOMEGA; j++) {
      std::cout << FGM[i][j] << "\t";
    }
    std::cout << std::endl;
  }
}
