#include <iostream>
#include <cmath>

#include <vector>
#include "rpoly.cpp"
#include "RunnerModule.hpp"
#include "simul_constants.hpp"

double max(double a, double b) {
  if (a > b) return a;
  else return b;
}

void coeff(double Omega, double dlnOmegadlnr, double kR, double kZ,
	   double ans[]){

  double dOmegadr, k2, k4, k6, k8, k10, A, B;
  double kva2, kva4;

  // recover dOmegadr w/ dlnOmegadlnr = r/Omega*dOmegadr
  dOmegadr = dlnOmegadlnr * Omega/simul::r;
  // precompute powers of k
  k2 = pow(kR,2) + pow(kZ,2);
  k4 = pow(k2,2);
  k6 = pow(k2,3);
  k8 = pow(k2,4);
  k10 = pow(k2,5);

  // Some elements :

  // D(u) = kR/kZ du/dZ - du/dR
  //      = (kR/kZ cos(O) - sin(O)) * du/dr + (...) * du/dO
  // the right term is null for Omega (and we assume for P also)

  // u_r     = sin theta u_R + cos theta u_Z
  // u_theta = cos theta u_R - sin theta u_Z

  // precompute 1/gamma*rho * ...
  A = ( 1/(simul::gamma*simul::rho) * pow((kR/kZ/cos(simul::theta) - 1/sin(simul::theta)),2) * (-simul::N2) );

  // precompute 
  //  1/pow(R,3) D(pow(R,4)*pow(Omega,2)) = 2*simul::*Omega*(kR/kZ - 1)*dOmega/dr-4*pow(Omega,2)
  B = ( 2 *simul::R*Omega* (kR/kZ/cos(simul::theta) - 1/sin(simul::theta)) * dOmegadr 
	- 4*pow(Omega,2) );

  // precompute (k * v_a)pow(,2) and its powers
  // k v_a = ksimul::theta * Bsimul::theta / 4pirho
  //       = Bsimul::theta / ... * (kR cos simul::theta - kZ sin simul::theta)
  kva2 = pow(simul::v_a_theta * (kR*cos(simul::theta) - kZ*sin(simul::theta)), 2);
  kva4 = pow(kva2,2);

  // compute the coefficients
  ans[0] = k2/pow(kZ,2);
  ans[1] = ans[0]*(2*simul::nu + 2*simul::eta + simul::xi)*k2;
  ans[2] = ( ans[0]*( k4*( pow(simul::nu,2) + pow(simul::eta,2) + 4*simul::nu*simul::eta + 
		   2*simul::nu*simul::xi + 2*simul::eta*simul::xi )
	      + 2*kva2)
	 - A - B );

  ans[3] = ( ans[0]*( ( 2*simul::eta*pow(simul::nu,2) + 2*simul::nu*pow(simul::eta,2) + pow(simul::nu,2)*simul::xi + pow(simul::eta,2)*simul::xi
		+ 4*simul::nu*simul::eta*simul::xi ) * k6
	      + 2*(simul::nu + simul::eta + simul::xi)*k2 * kva2 )
	 - (2*simul::eta + simul::nu)*k2 * A
	 - (2*simul::eta + simul::xi)*k2 * B );

  ans[4] = ( ans[0]*( (2*simul::eta*simul::xi*pow(simul::nu,2) + 2*simul::nu*pow(simul::eta,2)*simul::xi + pow(simul::eta,2)*pow(simul::nu,2))*k8
	      + 2*(simul::nu*simul::eta + simul::nu*simul::xi + simul::eta*simul::xi)*k4*kva2
	      + kva4)
	 - ( (2*simul::nu*simul::eta*k4 + pow(simul::eta,2)*k4 + kva2) * A )
	 - ( (2*simul::eta*simul::xi*k4 + pow(simul::eta,2)*k4 + kva2) * B )
	 - 4*pow(Omega,2)*kva2 );

  ans[5] = ( ans[0]*( simul::xi*pow(simul::eta,2)*pow(simul::nu,2)*k10 + 2*simul::xi*simul::nu*simul::eta*k6*kva2
	      + simul::xi*k2*kva4 )
	 - ( (simul::nu*pow(simul::eta,2)*k6 + simul::eta*k2*kva2) * A )
	 - ( (simul::xi*pow(simul::eta,2)*k6 + simul::xi*k2*kva2 ) * B )
	 - 4*pow(Omega,2)*kva2*simul::xi*k2 );
}
 


// Loop between Omega_range[om_b] & Omega_range[om_t] 
// (resp dom_b & dom_t for dlnOmegadlnr)
void loop(double FGMs[NOMEGA][NDOMEGA], double Omega_range[], double dlnOmegadlnr_range[],
	  double k_range[], int om_b, int om_e,
	  int dom_b, int dom_e, int nk) {
  
  double Omega, dlnOmegadlnr;
  double kR, kZ;
  double ret[3];
  // Iterate over all Omega, dlnOmega couples
  for (int a = om_b; a < om_e; a++) {
    for (int b = dom_b; b < dom_e; b++) {
      Omega = Omega_range[a];
      dlnOmegadlnr = dlnOmegadlnr_range[b];
      // max_kR = 0;
      // max_kZ = 0;

      get_FGM(Omega, dlnOmegadlnr, ret);
      FGMs[a][b] = ret[0];
      kR = ret[1];
      kZ = ret[2];
      

      if (simul::verbose) {
	std::cout << Omega << "\t" <<  dlnOmegadlnr 
		  << "\t" << FGMs[a][b] << "\t"
		  << "\t" << kR << "\t" << kZ << std::endl;
      }
    }
  }
}


void get_FGM(double Omega, double dlnOmegadlnr, double* ret) {
  double local_FGM = 0, FGM = 0;
  double kR, kZ;
  double kR_FGM_tmp, kZ_FGM_tmp;
  double polynomial[6];
  int order, nroots;
  double zeroi[5], zeror[5];
  int info[5];
  // iterate over all k-couples
  for (int c = 0; c < 2*NK; c++) {
    for (int d = 0; d < 2*NK; d++) {
      kR = simul::k_range[c];
      kZ = simul::k_range[d];

      coeff(Omega, dlnOmegadlnr, kR, kZ, polynomial);

      // Calculate the order of the polynomial
      // by looping from the end and decreasing as much as necessary
      order = 5;

      // Pass the polynomial to the rpoly function
      // result in zeror, zeroi
      nroots = rpoly(polynomial, order, zeror, zeroi, info);
            
      // local Fastest Growing Mode = biggest real value
      local_FGM = 0;
      for (int n = 0; n < nroots; n++) {
	local_FGM = max(local_FGM, zeror[n]);
      }
	  
      // if a bigger FGM, save infos
      if ( local_FGM > FGM ) {
	FGM = local_FGM;
	kR_FGM_tmp = kR;
	kZ_FGM_tmp = kZ;

	// Do some output if necessary
	if (simul::vverbose) {
	  std::cout << "# New FGM " << FGM << std::endl;
	}
      }
      if (simul::vverbose) {
	std::cout << Omega << "\t" << dlnOmegadlnr << "\t" <<
	  kR << "\t" <<  kZ << "\t" << local_FGM << std::endl;
      }
    }
  }
  ret[0] = FGM;
  ret[1] = kR_FGM_tmp;
  ret[2] = kZ_FGM_tmp;

}
