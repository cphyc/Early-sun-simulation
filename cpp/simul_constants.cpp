#include <cmath>
#include "simul_constants.hpp"

#ifdef USE_RG
// Red giant parameters
namespace simul {
  bool verbose = true, vverbose = false;

  // Constants
  double pi			= 3.1592653589793238462643383279502881971693;
  double Rsun			= 6.96e10;         //  cm
  double OmegaSun		= 2.7e-6;      //  s⁻¹
  double Rgas			= 8.314e7;         //  erg.K⁻¹.mol⁻¹
  double Msun			= 1.99e33;         //  g

  // set the parameters @ t	= 4.7e9yr:
  // Optional if using pressure scale
  // double g	= 1000;         //  cgs (taken at the surface)
  double nu		= 5.69;         //  cm².s⁻¹
  double eta		= 496;          //  cm².s⁻¹
  double xi		= 1.78e17;      //  cm².s⁻¹
  double r		= 2.07*Rsun;    //  cm
  double theta          = pi/4;         //  rda
  double R              = r*cos(theta); //  cm
  double gamma		= 1.63;         //  
  double rho		= 4.61e2;       //  g.cm⁻³
  double T		= 2.62e7;       //  K
  double B_theta	= 0;            //  Gauss (cgs)
  double v_a_theta      = B_theta/sqrt(4*pi*rho);
                                        // B_theta / np.sqrt(4*pi*rho)
                                        //  cm.s⁻¹
  double N2             = 2.22e-2;      //  s⁻²
  double H_pressure     = 3.79e-3*Rsun; //  cm

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
// Easly sun parameters
namespace simul {
  bool verbose = true, vverbose = false;

  double pi = 3.1592653589793238462643383279502881971693;
  double Rsun = 696342e5;;           //  cm
  double OmegaSun = 2.7e-6;          //  s⁻¹
  double Rgas = 8.314e7;             //  erg.K⁻¹.mol⁻¹
  double g = 27542.29;               //  cgs (taken at the surfce)

  // set the parameters :
  double nu = 15;                    //  cm².s⁻¹
  double eta = 496;                  //  cm².s⁻¹
  double xi = 4e6;                   //  cm².s⁻¹R
  double r = 0.7*Rsun;               //  cm
  double theta = pi/4;               //  rad
  double R = r*cos(theta);           //  cm
  double gamma = 5/3;                //  
  double rho = 0.4;                  //  g.cm⁻³
  double T = 2.6e6;                  //  K
  double B_theta = 0;                //  Gauss (cgs)
  double v_a_theta = B_theta/sqrt(4*pi*rho);
                                     // B_theta / np.sqrt(4*pi*rho)
                                     //  cm.s⁻¹
  double N2 = 6e-6;                  //  s⁻²

  // k ranges
  double lmax = 1e-2*Rsun;            //  cm
  double lmin = 1e-12*Rsun;           //  cm

  double k_range[2*NK];
  double Omega_range[NOMEGA];
  double dlnOmegadlnr_range[NDOMEGA];
}
#endif
