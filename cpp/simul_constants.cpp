#include <cmath>
#include "simul_constants.hpp"

namespace simul {
  bool verbose = true, vverbose = false;

  // Constants
  double pi = 3.1592653589793238462643383279502881971693;
  double Rsun = 6.96e10;         //  cm
  double OmegaSun = 2.7e-6;      //  s⁻¹
  double Rgas = 8.314e7;         //  erg.K⁻¹.mol⁻¹

  // set the parameters @ t = 4.7e9yr:
  // Optional if using pressure scale
  // double g = 1000;            //  cgs (taken at the surface)
  double nu = 5.71;              //  cm².s⁻¹
  double eta = 496;              //  cm².s⁻¹
  double xi = 2.13e16;           //  cm².s⁻¹
  double r = 2.97*Rsun;          //  cm
  double theta = pi/4;           //  rda
  double R = r*cos(theta);       //  cm
  double gamma = 1.63;           //  
  double rho = 4.59e2;           //  g.cm⁻³
  double T = 2.6e7;              //  K
  double B_theta = 4;            //  Gauss (cgs)
  double v_a_theta = B_theta/sqrt(4*pi*rho);
  // B_theta / np.sqrt(4*pi*rho)
                                 //  cm.s⁻¹
  double N2 = 4.59e-3;           //  s⁻²
  double H_pressure = 3.79e-3*Rsun;

  // k ranges
  // Either pressure scale height or Rgas*T/gx
  double lmax = H_pressure;      //  cm
  double lmin = 1e5;             //  cm

  double k_range[2*NK];
  double Omega_range[NOMEGA];
  double dlnOmegadlnr_range[NDOMEGA];
}
