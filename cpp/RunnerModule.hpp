#include "simul_constants.hpp"

// Compute a single FGM
double get_FGM(double, double);

// Loop between Omega_range[om_b] & Omega_range[om_t] 
// (resp dom_b & dom_t for dlnOmegadlnr)
void loop(double FGMs[NOMEGA][NDOMEGA], double Omega_range[],
	  double dlnOmegadlnr_range[],
	  double k_range[], int om_b, int om_e,
	  int dom_b, int dom_e, int nk);

void coeff(double Omega, double dlnOmegadlnr, double kR, double kZ,
	   double ans[]);

double get_FGM(double Omega, double dOmega);
