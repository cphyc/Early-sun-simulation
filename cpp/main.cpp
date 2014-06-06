#include <iostream>
#include "RunnerModule.hpp"
#include "simul_constants.hpp"
#include <cmath>

int main(){
  // To store the results
  double FGMs[NOMEGA][NDOMEGA];

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
  
  for (int i = 0; i < NOMEGA; i++) {
    for (int j = 0; j < NOMEGA; j++) {
      FGMs[i][j] = get_FGM(simul::Omega_range[i], simul::dlnOmegadlnr_range[j]);
      std::cout << FGMs[i][j] << "\t";
    }
    std::cout << "\n";
  }

}
