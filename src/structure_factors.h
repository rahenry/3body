#ifndef STRUCTURE_FACTORS_H
#define STRUCTURE_FACTORS_H

#include "general.h"

struct hyper_integrand_helper{
  vector<double> *transformed_state;
  double h;
  double eta;
  int kappa;
  int lr, nr, lrho, nrho;
  double rescale_param;

  hyper_integrand_helper(vector<double> *ts_, double h_, double eta_, int kappa_, int lr_, int nr_, int lrho_, int nrho_, double rescale_param_) : transformed_state(ts_), h(h_), eta(eta_), kappa(kappa_), lr(lr_), nr(nr_), lrho(lrho_), nrho(nrho_), rescale_param(rescale_param_){};
};

double hyper_helper_function(double rho, void * params);
double hyper_helper_function_scaled(double rho, void * params);

#endif

