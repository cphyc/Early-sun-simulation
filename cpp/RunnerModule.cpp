#include <iostream>
#include <cmath>

#include <vector>
#include "rpoly.cpp"
#include "RunnerModule.hpp"

bool verbose = true, vverbose = false;

Simul::Simul() {
  pi = 3.1592653589793238462643383279502881971693;
  Rsun = 696342e5;;           //  cm
  OmegaSun = 2.7e-6;          //  s⁻¹
  Rgas = 8.314e7;             //  erg.K⁻¹.mol⁻¹
  g = 27542.29;               //  cgs (taken at the surfce)

  // set the parameters :
  nu = 15;                    //  cm².s⁻¹
  eta = 0;                    //  cm².s⁻¹
  xi = 4e6;                   //  cm².s⁻¹
  r = 0.7*Rsun;               //  cm
  theta = pi/4;               //  rad
  R = r*cos(theta);           //  cm
  gamma = 5/3;                //  
  rho = 0.4;                  //  g.cm⁻³
  T = 2.6e6;                  //  K
  B_theta = 0;                //  Gauss (cgs)
  v_a_theta = 0;              // B_theta / np.sqrt(4*pi*rho)
  //  cm.s⁻¹
  N2 = 6e-6;                  //  s⁻²

  // k ranges
  lmax = Rgas*T/g;            //  cm
  lmin = 1e5;                 //  cm
}

double max(double a, double b) {
  if (a > b) return a;
  else return b;
}


void Simul::coeff(double Omega, double dlnOmegadlnr, double kR, double kZ,
	   double ans[]){

  double dOmegadr, k2, k4, k6, k8, k10, A, B;
  double kva2, kva4;

  // recover dOmegadr w/ dlnOmegadlnr = r/Omega*dOmegadr
  dOmegadr = dlnOmegadlnr * Omega/this->r;
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
  A = ( 1/(this->gamma*this->rho) * pow((kR/kZ/cos(this->theta) - 1/sin(this->theta)),2) * (-this->N2) );

  // precompute 
  //  1/pow(R,3) D(pow(R,4)*pow(Omega,2)) = 2*this->*Omega*(kR/kZ - 1)*dOmega/dr-4*pow(Omega,2)
  B = ( 2 *this->R*Omega* (kR/kZ/cos(this->theta) - 1/sin(this->theta)) * dOmegadr 
	- 4*pow(Omega,2) );

  // precompute (k * v_a)pow(,2) and its powers
  // k v_a = kthis->theta * Bthis->theta / 4pirho
  //       = Bthis->theta / ... * (kR cos this->theta - kZ sin this->theta)
  kva2 = this->v_a_theta * (kR*cos(this->theta) - kZ*sin(this->theta));
  kva4 = pow(kva2,2);

  // compute the coefficients
  ans[0] = k2/pow(kZ,2);
  ans[1] = ans[0]*(2*this->nu + 2*this->eta + this->xi)*k2;
  ans[2] = ( ans[0]*( k4*( pow(this->nu,2) + pow(this->eta,2) + 4*this->nu*this->eta + 
		   2*this->nu*this->xi + 2*this->eta*this->xi )
	      + 2*kva2)
	 - A - B );

  ans[3] = ( ans[0]*( ( 2*this->eta*pow(this->nu,2) + 2*this->nu*pow(this->eta,2) + pow(this->nu,2)*this->xi + pow(this->eta,2)*this->xi
		+ 4*this->nu*this->eta*this->xi ) * k6
	      + 2*(this->nu + this->eta + this->xi)*k2 * kva2 )
	 - (2*this->eta + this->nu)*k2 * A
	 - (2*this->eta + this->xi)*k2 * B );

  ans[4] = ( ans[0]*( (2*this->eta*this->xi*pow(this->nu,2) + 2*this->nu*pow(this->eta,2)*this->xi + pow(this->eta,2)*pow(this->nu,2))*k8
	      + 2*(this->nu*this->eta + this->nu*this->xi + this->eta*this->xi)*k4*kva2
	      + kva4)
	 - ( (2*this->nu*this->eta*k4 + pow(this->eta,2)*k4 + kva2) * A )
	 - ( (2*this->eta*this->xi*k4 + pow(this->eta,2)*k4 + kva2) * B )
	 - 4*pow(Omega,2)*kva2 );

  ans[5] = ( ans[0]*( this->xi*pow(this->eta,2)*pow(this->nu,2)*k10 + 2*this->xi*this->nu*this->eta*k6*kva2
	      + this->xi*k2*kva4 )
	 - ( (this->nu*pow(this->eta,2)*k6 + this->eta*k2*kva2) * A )
	 - ( (this->xi*pow(this->eta,2)*k6 + this->xi*k2*kva2 ) * B )
	 - 4*pow(Omega,2)*kva2*this->xi*k2 );
}



// Loop between Omega_range[om_b] & Omega_range[om_t] 
// (resp dom_b & dom_t for dlnOmegadlnr)
void Simul::loop(arr FGMs, double Omega_range[], double dlnOmegadlnr_range[],
		 double k_range[], int om_b, int om_e,
		 int dom_b, int dom_e, int nk) {
  
  double FGM, local_FGM, Omega, dlnOmegadlnr, kR, kZ;
  // double max_kR, max_kZ;
  double polynomial[6];
  int order, nroots;
  double zeroi[5], zeror[5];
  int info[5];

  // Iterate over all Omega, dlnOmega couples
  for (int a = om_b; a < om_e; a++) {
    for (int b = dom_b; b < dom_e; b++) {
      FGM = 0;
      Omega = Omega_range[a];
      dlnOmegadlnr = dlnOmegadlnr_range[b];
      // max_kR = 0;
      // max_kZ = 0;

      // iterate over all k-couples
      for (int c = 0; c < nk; c++) {
	for (int d = 0; d < nk; d++) {
	  kR = k_range[c];
	  kZ = k_range[d];

	  coeff(Omega ,dlnOmegadlnr, kR, kZ, polynomial);

	  // Calculate the order of the polynomial
	  // by looping from the end and decreasing as much as necessary
	  order = 5;
	  // for (int n = order; n >= 0; n--)  {
	  //   if (polynomial[n] == 0 ) 
	  //     order--;
	  //   else
	  //     break;
	  // }

	  // Pass the polynomial to the rpoly function
	  // result in zeror, zeroi
	  nroots = rpoly(polynomial, order, zeror, zeroi, info);
            
	  // local Fastest Growing Mode = biggest real value
	  local_FGM = 0;
	  for (int n = 0; n < nroots; n++) {
	    local_FGM = max(local_FGM, zeror[n]);
	  }
	  std::cout << kR << "\t" << kZ << "\t" << local_FGM << std::endl;
	  
	  // check whether we found a new absolute FGM
	  if ( local_FGM > FGM ) {
	    FGM = local_FGM;
	    // max_kR = kR;
	    // max_kZ = kZ;

	    // Do some output if necessary
	    if (vverbose) {
	      std::cout << "# New FGM " << FGM << std::endl;
	    }
	  }
	  if (vverbose) {
	    std::cout << Omega << "\t" << dlnOmegadlnr << "\t" <<
	      kR << "\t" <<  kZ << "\t" << local_FGM << std::endl;
	  }
	}
      }

      FGMs[a][b] = FGM/this->OmegaSun;

      if (verbose) {
	std::cout << Omega << "\t" <<  dlnOmegadlnr 
		  << "\t" << FGMs[a][b] << std::endl;
      }
    }
  }
}
