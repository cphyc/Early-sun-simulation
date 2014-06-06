#ifndef CONST_H
#define CONST_H

#include <cmath>
const int NOMEGA = 100, NDOMEGA = 100, NK = 100;

namespace simul {
  extern bool verbose, vverbose;

  extern double pi;
  extern double Rsun;
  extern double OmegaSun;
  extern double Rgas;
  extern double g;

  // set the parameters :
  extern double nu;
  extern double eta;
  extern double xi;
  extern double r;
  extern double theta;
  extern double R;
  extern double gamma;
  extern double rho;
  extern double T;
  extern double B_theta;
  extern double v_a_theta;
  //  cm.s⁻¹
  extern double N2;

  // k ranges
  extern double lmax;
  extern double lmin;

  extern double k_range[NK];
  extern double Omega_range[NOMEGA];
  extern double dlnOmegadlnr_range[NDOMEGA];
};
#endif
