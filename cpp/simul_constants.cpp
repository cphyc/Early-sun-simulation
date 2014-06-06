#include <cmath>
#include "simul_constants.hpp"

namespace simul {
  bool verbose = true, vverbose = false;

  double pi = 3.1592653589793238462643383279502881971693;
  double Rsun = 696342e5;;           //  cm
  double OmegaSun = 2.7e-6;          //  s⁻¹
  double Rgas = 8.314e7;             //  erg.K⁻¹.mol⁻¹
  double g = 27542.29;               //  cgs (taken at the surfce)

  // set the parameters :
  double nu = 15;                    //  cm².s⁻¹
  double eta = 0;                    //  cm².s⁻¹
  double xi = 4e6;                   //  cm².s⁻¹
  double r = 0.7*Rsun;               //  cm
  double theta = pi/4;               //  rad
  double R = r*cos(theta);           //  cm
  double gamma = 5/3;                //  
  double rho = 0.4;                  //  g.cm⁻³
  double T = 2.6e6;                  //  K
  double B_theta = 0;                //  Gauss (cgs)
  double v_a_theta = 0;              // B_theta / np.sqrt(4*pi*rho)
  //  cm.s⁻¹
  double N2 = 6e-6;                  //  s⁻²

  // k ranges
  double lmax = Rgas*T/g;            //  cm
  double lmin = 1e5;                 //  cm

  double k_range[NK];
  double Omega_range[NOMEGA];
  double dlnOmegadlnr_range[NDOMEGA];
}
