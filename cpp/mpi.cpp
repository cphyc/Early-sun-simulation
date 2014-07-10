#include "mpi.h"
#include <iostream>
#include "RunnerModule.hpp"
#include "simul_constants.hpp"
#define CONTINUE 1
#define STOP 0
#define MANAGER 0
// SHIFT is used to communicate between the worker and the manager
#define SHIFT 12

// verbose 1 : output everything, 2  %
int verbosity = 2;

// output arrays !
double FGM[NOMEGA][NDOMEGA];
double kR_FGM[NOMEGA][NDOMEGA];
double kZ_FGM[NOMEGA][NDOMEGA];

// First call to next_pos will set posdOmega->0
int posOmega = 0, posdOmega = -1; 


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
void print_array(double FGM[NOMEGA][NDOMEGA]);


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

/* Modifies the position to the next position in the list and return True
   or returns False if nothing possible anymore */
bool next_pos() {
  // Increment for dOmega ok ?
  if (++posdOmega < NDOMEGA)
    return true;
  else {
    posdOmega = 0;
    if (++posOmega < NOMEGA)
      return true;
    else
      return false;
  }
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
  int om_index_ret, sender; 
  int dom_index_ret;
  double ret[3];
  int ins[2];
  int i;
  
  // number of jobs to launch
  int njobs = min(numprocs - 1, NOMEGA*NDOMEGA);

  if (verbosity == 1)// Send a first row
    printf("M : %d job%s\n", njobs, (njobs > 1)? "s" : "");
  
  i=1;
  while (i <= njobs && next_pos()) {
    if (verbosity == 1) {
      printf("M : sending first instructions to %d with %d-%d...\n", i,
	     posOmega, posdOmega);
    }
    // Instructions = indexes of Omega, dOmega
    ins[0] = posOmega;
    ins[1] = posdOmega;

    MPI_Send(ins,               // input
	     2,                 // size
	     MPI_INT,           // type
	     i,                 // worker #i
	     CONTINUE,          // continue
	     MPI_COMM_WORLD);   // communicator
    i++;
  }
  
  /* For all the elements in the matrix, get an answer and send a next one
     if available */
  for (int i = 0; i < NOMEGA*NDOMEGA; i++) {
    if (verbosity == 1)
      printf("M : waiting for an answer ... ");

    MPI_Recv(ret,
	     3,
	     MPI_DOUBLE,
	     MPI_ANY_SOURCE,
	     MPI_ANY_TAG,
	     MPI_COMM_WORLD,
	     &status);

    /* Fetch back the omega rank and the sender.*/
    sender  = status.MPI_SOURCE;
    /* om_index is in the first SHIFT bits on the left */
    om_index_ret = status.MPI_TAG >> SHIFT;
    /* dom_index is on the first SHIFT bits on the right */
    dom_index_ret = status.MPI_TAG % ( 1 << SHIFT );
    
    if (verbosity == 1) {
      printf("got %d-%d from %d, storing it ... ", om_index_ret, dom_index_ret,
	     sender);
    }
    else if (verbosity == 2){ 
      printf("%.2f%%", (NOMEGA*posOmega + posdOmega) * 100. / (NOMEGA*NDOMEGA));
      std::cout << std::endl;
    }

    /* Store the answer in our array */
    FGM[om_index_ret][dom_index_ret] = ret[0];
    kR_FGM[om_index_ret][dom_index_ret] = ret[1];
    kZ_FGM[om_index_ret][dom_index_ret] = ret[2];

    if (verbosity == 1)
      printf("stored. \n");

    /* if next position is available, send it */
    if (next_pos()){
      if (verbosity == 1){
	printf("M : sending instructions to %d with %d-%d ...", sender,
	       posOmega, posdOmega);
      }

      ins[0] = posOmega;
      ins[1] = posdOmega;
      MPI_Send(ins,
	       2,
	       MPI_INT,
	       sender,
	       CONTINUE,
	       MPI_COMM_WORLD);
      if (verbosity == 1)
	printf(" sent!\n");
    }
    /* No other works, tell the sender to stop by sending the tag -1*/
    else {
      if (verbosity == 1)
	printf("M : No more work for %d :(.\n", sender);
      MPI_Send(MPI_BOTTOM,
	       0,
	       MPI_DOUBLE,
	       sender,
	       STOP,
	       MPI_COMM_WORLD);
    }
  }
  if (verbosity == 1)
    printf("M : Printing out the total answer\n");

  print_array(FGM);
  print_array(kR_FGM);
  print_array(kZ_FGM);
  return;
}


void worker_code(void) {
  MPI_Status status;
  int rank;
  double Omega, dOmega;
  double FGM;
  int om_index, dom_index;
  int index[2];
  unsigned int flag;
  double ret[3];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  /* Get a new command from 0 with a null (0) tag. The size of the
     input in NDOMEGA, the last one being Omega*/

  MPI_Recv(index,
	   2,
	   MPI_INT,
	   MANAGER,
	   MPI_ANY_TAG, 
	   MPI_COMM_WORLD, &status);
    
  while (status.MPI_TAG == CONTINUE) {
    om_index = index[0];
    dom_index = index[1];
      
    /* Call the get_FGM function that finds the fastest 
       growing mode at a given Omega by looping over dOmega */
    Omega  = simul::Omega_range[om_index];
    dOmega = simul::dlnOmegadlnr_range[dom_index];
    if (verbosity == 1)
      printf("W%d: working on %d - %d (%e - %e)\n", rank, om_index, dom_index, Omega, dOmega);
    
    get_FGM(Omega, dOmega, ret);

    /* l-t-r :  SHIFT bits for om_index _  SHIFT bits for dom_index
       ex : nOmega = 5, ndOmega = 10, SHIFT = 8
       
       flag = 00000101 00001010
            =   5*2^8 +      10   */
    flag = (om_index << SHIFT) | dom_index;

    /* Send the answer (42) to 0 with the Omega_rank.*/
    MPI_Send(ret,
	     3,
	     MPI_DOUBLE,
	     MANAGER,
	     flag,
	     MPI_COMM_WORLD);

    /* Get a new job */
    MPI_Recv(index,
	     2,
	     MPI_INT,
	     MANAGER,
	     MPI_ANY_TAG, 
	     MPI_COMM_WORLD, &status);
  }
  printf("W%d: My work is over, I shall rest in peace.\n", rank);
  return ;
}

void print_array(double FGM[NOMEGA][NDOMEGA]){
  for (int i = 0; i < NOMEGA; i++) {
    for (int j = 0; j < NDOMEGA; j++) {
      std::cout << FGM[i][j] << "\t";
    }
    std::cout << std::endl;
  }
}
