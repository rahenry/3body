#include "general.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include "integration.h"
#include "findV0.h"

double calc_V0(double r0, double a0){
  V0Finder_new finder;
  finder.findV0(r0, 1.0/a0);
  return finder.V0;
}

double root_function(double x, void *params){
  V0Finder_new *V0F = (V0Finder_new *) params;
  V0F->V0 = x;
  V0F->finda();
  return 1.0 / V0F->a - V0F->ainvDesired;

}

void V0Finder_new::findV0(double r0_, double ainvDesired_){
  ainvDesired = ainvDesired_;
  gsl_set_error_handler_off();
  int status, iter=0, max_iter=100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  gsl_function F;
  F.function = &root_function;
  F.params = this;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);

  r0 = r0_;
  V0Increment = 30.0;

  if (r0 > 0.01) V0 = -300.0;
  else V0 = -3000.0;
  finda();
  double a0 = a;
  int sign0 = (1.0/a0 - ainvDesired > 0) - (1.0/a0 - ainvDesired < 0);

  ofstream errorfile("V0errorfile");
  errorfile.precision(16);
  errorfile << sign0 << endl;

  bool check = false;
  while (!check){
    finda();
    errorfile << V0 << " " << findainv() << endl;
    if ( ((1.0/a - ainvDesired > 0) - (1.0/a - ainvDesired < 0)) != sign0 ){
      check = true;
    }
    else{
      V0Increment *= 1.05;
      V0 -= V0Increment;
    }
  };

  //cout << "V0 test: " << V0 << " " << V0Increment << endl;
  double V0_low = V0 - V0Increment;
  double V0_high = V0 + V0Increment;

  status = gsl_root_fsolver_set(s, &F, V0_low, V0_high);

  do{
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r = gsl_root_fsolver_root(s);
    V0_low = gsl_root_fsolver_x_lower(s);
    V0_high = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(V0_low, V0_high, 0, 1E-8);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free(s);
  V0 = r;
}

void V0Finder_new::findV0(double r0_){
  findV0(r0_, 0.0);
}

int integration_helper_u(double x, const double u[], double f[], void *params){
  V0Finder_new *V0F = (V0Finder_new *) params;
  f[0] = u[1];
  f[1] = u[0] * V0F->V0 * exp(-x * x / V0F->r0 / V0F->r0 / 2.0 - V0F->k * V0F->k);
  //cout << u[0] << " " << u[1] << " " << f[0] << " " << f[1] << " | " << V0F->r0 << " " << V0F->k << " " << V0F->V0 << endl;
  return GSL_SUCCESS;
}

int jacob(double x, const double u_old[], double *dfdu, double dfdx[], void *params){
  gsl_matrix_view dfdu_mat = gsl_matrix_view_array(dfdu, 2, 2);
  gsl_matrix *m = &dfdu_mat.matrix;
  V0Finder_new *V0F = (V0Finder_new *) params;
  gsl_matrix_set(m, 0, 0, 1.0);
  gsl_matrix_set(m, 0, 1, 1.0);
  gsl_matrix_set(m, 1, 1, 0.0);
  gsl_matrix_set(m, 1, 0, V0F->V0 * exp(-x * x / V0F->r0 / V0F->r0 - V0F->k * V0F->k));
  dfdx[0] = 0.0;
  dfdx[1] = u_old[0] * V0F->V0 * -1.0 * x / V0F->r0 / V0F->r0 * exp(-x * x / V0F->r0 / V0F->r0 - V0F->k * V0F->k);

  cout << dfdx[0] << " " << dfdx[1] << " " << dfdu[0] << endl;
  return GSL_SUCCESS;
}

void V0Finder_new::findu(){
  u[0] = 0.0;
  u[1] = 1.0;

  gsl_odeiv2_system sys = {integration_helper_u, jacob, 2, this};
  gsl_odeiv2_driver *d = 
    gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1E-6, 1E-6, 0.0);
  double x = 0.0;

  for (int i = 1; i<=5; i++){
    double xi = i * rMax / 5.0;
    int status = gsl_odeiv2_driver_apply(d, &x, xi, &u[0]);

    if (status != GSL_SUCCESS){
      cout << "NO " << endl;
      break;
    }
  }

  gsl_odeiv2_driver_free(d);
};

void V0Finder_new::finda(){
  double aMax = 1.0, aMaxPrev = 0.0;
  int index=1;
  while (abs(aMax-aMaxPrev) > 1E-8){
    index++;
    k = pow(10.0,-index);
    init_wavenumber();
    findu();
    beta = rMax * u[1] / u[0] - 1.0;
    aMaxPrev = aMax;
    aMax = atan2(rMax * j0P - beta * j0, rMax * n0P - beta * n0) / k;
    if (index > 15) break;
  }
  a = aMax;
};

double V0Finder_new::findainv(){
  finda();
  return 1.0 / a;
}

void V0Finder_new::operator() (const vector<double> &u_current, vector<double> &dudx, const double x){
  dudx[0] = u_current[1];
  dudx[1] = (V0 * exp(-x*x/r0/r0/2.0) - 1*k*k) * u_current[0];
};


void V0Finder_new::init_wavenumber(){
  kR = k * rMax;   
  j0 = sin(kR) / k;
  j0P = cos(kR) - sin(kR) / kR; 
  n0 = -cos(kR) / k;  
  n0P = sin(kR) + cos(kR) / kR;
}

bool tol_1(double min, double max){
  bool res;
  //if ( 1.0 * abs(min-max) / (abs(min) + abs(max)) < 1E-16 
      //&& abs(min-max) < 1E-16){
  if ( 1.0 * abs(abs(min)-abs(max)) / (abs(min) + abs(max)) < 1E-16 
      && abs(abs(min)-abs(max)) < 1E-16){
    res = true;
  }
  else{
    res = false;
  }
  return res;
}
