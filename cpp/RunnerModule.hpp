#include <array>
const int N = 100;
typedef std::array< std::array<double, N>, N > arr;

class Simul {
public:
  int nk, nOmega, ndOmega;
  double pi, Rsun, OmegaSun, Rgas, g, nu, eta, xi;
  double r, theta, R, gamma, rho, T, B_theta, v_a_theta;
  double N2, lmax, lmin;

  // Set the parameters to default values
  Simul();

  // Loop between Omega_range[om_b] & Omega_range[om_t] 
  // (resp dom_b & dom_t for dlnOmegadlnr)
  void loop(arr FGMs, double Omega_range[], double dlnOmegadlnr_range[],
	    double k_range[], int om_b, int om_e,
	    int dom_b, int dom_e, int nk);

  void coeff(double Omega, double dlnOmegadlnr, double kR, double kZ,
	     double ans[]);
  
};
