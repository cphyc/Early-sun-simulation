#include <cmath>
#include "simul_constants.hpp"
#define USE_RG
#ifdef USE_RG
// Red giant parameters
namespace simul {
  //nu = 3.52e+02 	 eta = 6.51e+02  xi = 2.90e+16 	 r = 7.47e-01 	 gamma = 1.67e+00
  //rho = 9.01e-02 	 T = 3.29e+06 	 N^2_mean = 1.58e-07   	 H = 4.80e-01
  bool verbose = true, vverbose = false;

  // Constants
  double pi			= 3.1592653589793238462643383279502881971693;
  double Rsun			= 6.96e10;         //  cm
  double OmegaSun		= 2.7e-6;          //  s⁻¹
  double Rgas			= 8.314e7;         //  erg.K⁻¹.mol⁻¹
  double Msun			= 1.99e33;         //  g 
  double theta                  = pi/4;
  
  double nu		 =   2.69e-14; // 3.52e2 cms
  double eta		 =   4.98e-14; // 6.51e2 cms
  double xi		 =   2.22e0 ; // 2.90e16 cms
  double r		 =   7.47e-1;
  double gamma		 =   1.67e0;
  double rho		 =   9.01e2;
  double T		 =   3.29e6;
  double N2		 =   2.15e4; // 1.57e-7 cms
  double H_pressure	 =  4.80e-1;

  double B_theta	= 0;            //  Gauss (cgs)
  double v_a_theta      = B_theta/sqrt(4*pi*rho);
                                        // B_theta / np.sqrt(4*pi*rho)

  // k ranges
  //double lmax = Rgas*T/g;            //  cm
  //double lmin = 1e5;                 //  cm
  double lmax = 1e-2*Rsun;
  double lmin = 1e-14*Rsun;

  double k_range[2*NK];
  double Omega_range[NOMEGA];
  double dlnOmegadlnr_range[NDOMEGA];
}
#else
// Early sun parameters (non dimensionals !)
namespace simul {
  bool verbose = true, vverbose = false;

  double pi = 3.1415926536;
  double Rsun = 6.96e10;;        //  cm
  double OmegaSun = 2.7e-6;      //  cm^-1
  // double Rgas = 8.314e7;      //  erg.K⁻¹.mol⁻¹
  // double g = 27542.29;        //  cgs (taken at the surfce)

  // set the parameters :
  double nu = 1.15e-15;          //  Rsun**2*OmegaSun
  double eta = 3.79e-14;         //  Rsun**2*OmegaSun
  double xi = 3.06e-10;          //  Rsun**2*OmegaSun
  double r = 0.7;                //  Rsun
  double theta = pi/4.;          //  rad
  // double R = r*cos(theta);    //  cm
  // double gamma = 5/3;         //  
  double rho = 0.4;              //  g.cm⁻³
  // double T = 2.6e6;           //  K
  double B_theta = 0;            //  Gauss (cgs)
  double v_a_theta = B_theta/sqrt(4*pi*rho);
                                 // B_theta / np.sqrt(4*pi*rho)
  double N2 = 1.78e5;            //  in OmegaSun**2 unit

  // k ranges
  double lmax = 1e-2;            //  Rsun
  double lmin = 1e-14;           //  Rsun

  double k_range[2*NK];
  double Omega_range[NOMEGA];
  double dlnOmegadlnr_range[NDOMEGA];
}
#endif
