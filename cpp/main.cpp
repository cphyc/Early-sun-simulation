#include <iostream>
#include "RunnerModule.hpp"
#include <cmath>

int main(){
  // To store the results
  arr FGMs;
  
  // To create all the parameters and get the methods (see RunnerModule.{hpp,cpp})
  Simul s;

  // resolution
  int nk = 5;
  int nOmega = 5;
  int ndOmega = nOmega;
  double Omega_range[1000], dlnOmegadlnr_range[1000], k_range[1000];
  
  // logarithmic scale for k
  double kmin = 2*s.pi/s.lmax;
  double kmax = 2*s.pi/s.lmin;
  double fact = pow(kmax/kmin, 1./(nk-1));
  for (int n = 0; n < nk; n++) {
    double tmp = kmin * pow(fact, n);
    k_range[n]    = tmp;
    //    k_range[nk+n] = -tmp;
  }

  // Omega range and dlnOmegadlnr_range
  for (int n = 0; n < nOmega; n++)
    Omega_range[n] = 31./nOmega*s.OmegaSun*n;
  for (int n = 0; n < ndOmega; n++)
    dlnOmegadlnr_range[n] = -2.5/ndOmega*n;
  
  // // Put only zeros in FGMs
  // for (int n1 = 0; n1 < N; n1++){
  //   for (int n2 = 0; n2 < N; n2++){
  //     FGMs[n1][n2] = 0;
  //   }
  // }
  
  // Actually do the simulation
  s.loop(FGMs, Omega_range, dlnOmegadlnr_range, k_range,
	 0, nOmega,
	 0, ndOmega,
	 nk);

  // Print the output
  std::cout << "# Printing out FGMs" << std::endl;
  for (int n1 = 0; n1 < N; n1++){
    for (int n2 = 0; n2 < N; n2++){
      std::cout << FGMs[n1][n2] << "\t";
    }
    std::cout << std::endl;
  }
}
