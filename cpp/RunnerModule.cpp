#include <iostream>
#include <cmath>

#include <vector>
#include "rpoly.cpp"
#include "RunnerModule.hpp"
#include "simul_constants.hpp"
bool verbose_coeffs = false;
double max(double a, double b) {
  if (a > b) return a;
  else return b;
}

void coeff(double Omega, double dlnOmegadlnr, double kr, double kthe,
	   double ans[]){
  double k2, k4, k6, k8, k10, strat_term, rot_term;
  double kva2, kva4;
  double kR, kZ;

  kR = sin(simul::theta)*kr + cos(simul::theta)*kthe/simul::r;
  kZ = cos(simul::theta)*kr - sin(simul::theta)*kthe/simul::r;

  // recover dOmegadr w/ dlnOmegadlnr = r/Omega*dOmegadr
  // dOmegadr = dlnOmegadlnr * Omega/simul::r;
  // precompute powers of k
  k2 = pow(kr,2) + pow(kthe,2);
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


  double sthe = sin(simul::theta);
  double cthe = cos(simul::theta);

  // precompute 1/gamma*rho * ...
  strat_term = pow((kR/kZ*cthe - sthe), 2) * simul::N2;

  // precompute 
  //  1/pow(R,3) D(pow(R,4)*pow(Omega,2)) = 2*simul::*Omega*(kR/kZ - 1)*dOmega/dr-4*pow(Omega,2)
  // rot_term = ( 2 *simul::R*Omega* (kR/kZ*cthe - sthe) * dOmegadr 
  // rot_term = pow(Omega, 2) * 
  //   ( 2 * sthe * (kR/kZ*cthe - sthe) * dlnOmegadlnr 
  //     - 4 );
  rot_term=-1.*( (kR/kZ*cthe*sthe - sthe*sthe)*
		 (4.*Omega*Omega + 2*Omega*Omega*dlnOmegadlnr) 
		 - (kR/kZ*sthe+cthe)*(4.*Omega*Omega*cthe) );
  // R = r*sin(theta)
  // Om/r * (dlog Om / dlog r) = (d Om / d r)

  // precompute (k * v_a)pow(,2) and its powers
  // k v_a = ktheta * Btheta / 4pirho
  //       = Btheta / ... * (kR cos simul::theta - kZ sin simul::theta)
  kva2 = pow(simul::v_a_theta * kthe, 2);
  kva4 = pow(kva2, 2);

  // compute the coefficients
  ans[0] = k2/(kZ*kZ);
  ans[1] = ans[0]*(2*simul::nu + 2*simul::eta + simul::xi)*k2;
  ans[2] = ( ans[0]*( k4*( pow(simul::nu,2) 
			   + pow(simul::eta,2) 
			   + 4*simul::nu*simul::eta 
			   + 2*simul::nu*simul::xi 
			   + 2*simul::eta*simul::xi )
		      + 2*kva2 )
	     + strat_term + rot_term );

  ans[3] = ( ans[0]*( ( 2*simul::eta*pow(simul::nu,2) 
			+ 2*simul::nu*pow(simul::eta,2) 
			+ pow(simul::nu,2)*simul::xi 
			+ pow(simul::eta,2)*simul::xi
			+ 4*simul::nu*simul::eta*simul::xi ) * k6
		      + 2*(simul::nu + simul::eta + simul::xi)*k2 * kva2 )
	     + (2*simul::eta + simul::nu)*k2 * strat_term
	     + (2*simul::eta + simul::xi)*k2 * rot_term );

  ans[4] = ( ans[0]*( ( 2*simul::eta*simul::xi*pow(simul::nu,2) 
			+ 2*simul::nu*pow(simul::eta,2)*simul::xi 
			+ pow(simul::eta,2)*pow(simul::nu,2))*k8
		      + 2*( simul::nu*simul::eta 
			    + simul::nu*simul::xi 
			    + simul::eta*simul::xi)*k4*kva2
		      + kva4 )
	     + ( (2*simul::nu*simul::eta*k4 
		  + pow(simul::eta,2)*k4 
		  + kva2) * strat_term )
	     +( (2*simul::eta*simul::xi*k4 
		 + pow(simul::eta,2)*k4 
		 + kva2) * rot_term )
	     - 4*pow(Omega,2)*kva2 );

  ans[5] = ( ans[0]*( simul::xi*pow(simul::eta,2)*pow(simul::nu,2)*k10 
		      + 2*simul::xi*simul::nu*simul::eta*k6*kva2
		      + simul::xi*k2*kva4 )
	     + ( (simul::nu*pow(simul::eta,2)*k6 
		  + simul::eta*k2*kva2) * strat_term )
	     + ( (simul::xi*pow(simul::eta,2)*k6 
		  + simul::xi*k2*kva2 ) * rot_term )
	     - 4*pow(Omega,2)*kva2*simul::xi*k2 );

  if (verbose_coeffs) {
    printf("%2.5e\t%2.5e\t", kr, kthe);
    for (int i = 5; i >=0 ; i --)
      printf("%2.5e\t", ans[i]);
    std::cout << std::endl;
  }
}
 


// Loop between Omega_range[om_b] & Omega_range[om_t] 
// (resp dom_b & dom_t for dlnOmegadlnr)
void loop(double FGMs[NOMEGA][NDOMEGA], double Omega_range[], double dlnOmegadlnr_range[],
	  double k_range[], int om_b, int om_e,
	  int dom_b, int dom_e, int nk) {
  
  double Omega, dlnOmegadlnr;
  double kr, kthe;
  double ret[3];
  // Iterate over all Omega, dlnOmega couples
  for (int a = om_b; a < om_e; a++) {
    for (int b = dom_b; b < dom_e; b++) {
      Omega        = Omega_range[a];
      dlnOmegadlnr = dlnOmegadlnr_range[b];

      get_FGM(Omega, dlnOmegadlnr, ret);
      FGMs[a][b] = ret[0];
      kr         = ret[1];
      kthe       = ret[2];
      

      if (simul::verbose) {
	std::cout << Omega << "\t" <<  dlnOmegadlnr 
		  << "\t" << FGMs[a][b] << "\t"
		  << "\t" << kr << "\t" << kthe << std::endl;
      }
    }
  }
}


void get_FGM(double Omega, double dlnOmegadlnr, double* ret) {
  double local_FGM = 0, FGM = 0;
  double kr_FGM_tmp, kthe_FGM_tmp;
  double kr, kthe;
  double polynomial[6];
  int order, nroots;
  double zeroi[5], zeror[5];
  int info[5];
  // iterate over all k-couples
  for (int c = 0; c < 2*NK; c++) {
    for (int d = 0; d < 2*NK; d++) {
      kr   = simul::k_range[c];
      kthe = simul::k_range[d];

      coeff(Omega, dlnOmegadlnr, kr, kthe, polynomial);

      // Calculate the order of the polynomial
      // by looping from the end and decreasing as much as necessary
      order = 5;
      // if (polynomial[5] == 0) {
      // 	order --;
      // 	if (polynomial[4] == 0) {
      // 	    order --;
      // 	  }
      // }

      // Pass the polynomial to the rpoly function
      // result in zeror, zeroi
      nroots = rpoly(polynomial, order, zeror, zeroi, info);
            
      // local Fastest Growing Mode = biggest real value
      local_FGM = 0;
      for (int n = 0; n < 5; n++) {
	local_FGM = max(local_FGM, zeror[n]);
      }
	  
      // if a bigger FGM, save infos
      if ( local_FGM > FGM ) {
	FGM = local_FGM;
	// Please note that if FGM == 0, kR & kZ have no meaning.
	kr_FGM_tmp   = kr;
	kthe_FGM_tmp = kthe;

	// Do some output if necessary
	if (simul::vverbose) {
	  std::cout << "# New FGM " << FGM << std::endl;
	}
      }
      if (simul::vverbose) {
	std::cout << Omega << "\t" << dlnOmegadlnr << "\t" <<
	  kr << "\t" <<  kthe << "\t" << local_FGM << std::endl;
      }
    }
  }
  ret[0] = FGM;
  ret[1] = kr_FGM_tmp;
  ret[2] = kthe_FGM_tmp;
  return;
}
