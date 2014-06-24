#include "mpi.h"
#include <iostream>
#include "RunnerModule.hpp"
#include "simul_constants.hpp"
#define CONTINUE 1
#define STOP 0
#define MANAGER 0

// output array !
double FGM[NOMEGA][NDOMEGA];


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
    simul::dlnOmegadlnr_range[n] = -5./NDOMEGA*n;
}

int main(int argc, char *argv[]) {
  int numprocs, rank;

  // Disable buffering for stdout
  // setbuf(stdout, NULL);

  set_ranges();

  MPI_Init(&argc, &argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if (rank == MANAGER)
    manager_code(numprocs);
  else
    worker_code();

  MPI_Finalize();
  return 0;
}

void manager_code(int numprocs){
  MPI_Status status;
  int om_index = 0;    // This is the index of the max Omega that is calculated
  int om_index_ret, sender; // This is th
  double ret[NDOMEGA];
  
  // number of jobs to launch
  int njobs = min(numprocs - 1, NOMEGA);

  // Send a first row
  printf("M : %d job%s\n", njobs, (njobs > 1)? "s" : "");
  for (int i = 1 ; i <= njobs ; i++) {
    printf("M : sending first instructions to %d with omega index %d...", i, om_index);
    MPI_Send(&om_index,                   // input
	     1,                           // size
	     MPI_INT,                     // type
	     i,                           // worker #
	     CONTINUE,                    // continue
	     MPI_COMM_WORLD);             // communicator
    printf(" sent!\n");
    om_index++;
  }

  for (int i = 0; i < NOMEGA; i++) {
    printf("M : waiting for an answer ... ");

    MPI_Recv(ret,
	     NDOMEGA,
	     MPI_DOUBLE,
	     MPI_ANY_SOURCE,
	     MPI_ANY_TAG,
	     MPI_COMM_WORLD,
	     &status);

    /* Fetch back the omega rank and the sender.*/
    sender  = status.MPI_SOURCE;
    om_index_ret = status.MPI_TAG;

    printf("got column %d from %d, storing it ... ", om_index_ret, sender);

    /* Store the answer in our array */
    for (int dom_index = 0; dom_index < NDOMEGA; dom_index++) {
      FGM[om_index_ret][dom_index] = ret[dom_index];
    }
    printf("stored. \n");

    if (om_index < NOMEGA){
      printf("M : sending instructions to %d with omega index %d ...", sender, om_index);
      MPI_Send(&om_index,
	       1,
	       MPI_INT,
	       sender,
	       CONTINUE,
	       MPI_COMM_WORLD);
      printf(" sent!\n");
      om_index ++;
    }
    /* No other works, tell the sender to stop by sending the tag -1*/
    else {
      printf("M : No more work for %d :(.\n", sender);
      MPI_Send(MPI_BOTTOM,
	       0,
	       MPI_DOUBLE,
	       sender,
	       STOP,
	       MPI_COMM_WORLD);
    }
  }
  printf("M : Printing out the total answer\n");
  print_FGM(FGM);
  return;
}


void worker_code(void) {
  MPI_Status status;
  int rank;
  double Omega, dOmega;
  double FGM_col[NDOMEGA];
  int om_index;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  /* Get a new command from 0 with a null (0) tag. The size of the
     input in NDOMEGA, the last one being Omega*/

  MPI_Recv(&om_index,
	   1,
	   MPI_INT,
	   MANAGER,
	   MPI_ANY_TAG, 
	   MPI_COMM_WORLD, &status);
    
  while (status.MPI_TAG == CONTINUE) {
    printf("W%d: working on column %d\n", rank, om_index);
      
    /* Call the get_FGM function that finds the fastest 
       growing mode at a given Omega by looping over dOmega */
    for (int i = 0; i < NDOMEGA; i++) {
      Omega  = simul::Omega_range[om_index];
      dOmega = simul::dlnOmegadlnr_range[i];

      FGM_col[i] = get_FGM(Omega, dOmega);

      printf("W%d: working on %e, %e\n", rank, Omega, dOmega);
    }
      
    /* Send the answer (42) to 0 with the Omega_rank.*/
    MPI_Send(FGM_col,
	     NDOMEGA,
	     MPI_DOUBLE,
	     MANAGER,
	     om_index,
	     MPI_COMM_WORLD);

    /* Get a new job */
    MPI_Recv(&om_index,
	     1,
	     MPI_INT,
	     MANAGER,
	     MPI_ANY_TAG, 
	     MPI_COMM_WORLD, &status);
  }
  printf("W%d: My work is over, I shall rest in peace.\n", rank);
  return ;
}

void print_FGM(double FGM[NOMEGA][NDOMEGA]){
  for (int i = 0; i < NOMEGA; i++) {
    for (int j = 0; j < NDOMEGA; j++) {
      std::cout << FGM[i][j] << "\t";
    }
    std::cout << std::endl;
  }
}
