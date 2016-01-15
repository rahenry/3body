#ifndef FIND_V0
#define FIND_V0

#include "general.h"

double calc_V0(double r0, double a0);
double root_function(double x, void * params);
int integration_helper_u(double x, const double u[], double f[], void *params);
int jacob(double x, const double u_old[], double *dfdu, double dfdx[], void *params);


class V0Finder_new{
public:
  double r0, V0, rMin, rMax, k, beta, a, kR, j0, j0P, n0, n0P, V0Increment, ainvDesired;
  size_t steps;
  pair<double, double> result;
  vector<double> u;
  //typedef runge_kutta_cash_karp54< vector<double> > error_stepper_type;
  V0Finder_new() : rMin(0.0), rMax(1000.0), k(1E-10) {u.resize(2); kR = k * rMax; j0 = sin(kR) / kR; j0P = cos(kR) / kR - sin(kR) / kR / kR; n0 = -cos(kR) / kR; n0P = sin(kR) / kR + cos(kR) / kR / kR;};

  void findV0(double r0_);
  void findV0(double r0_, double ainvDesired_);
  void init_wavenumber();
  void findu();
  void finda();
  double findainv();
    
  void operator() (const vector<double> &u_current, vector<double> &dudx, const double x);
     
};
  
bool tol_1(double min, double max);
#endif
