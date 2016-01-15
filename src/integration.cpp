#include "integration.h"

double integrate_to_inf(gsl_function &F, double lower, double error_abs, double error_rel){
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000);
  double result, error;
  gsl_integration_qagiu(&F, lower, error_abs, error_rel, 100000, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

double integrate(gsl_function &F, double lower, double upper, double error_abs, double error_rel){
  int s = 10000;
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (s);

  double result, error;
  //gsl_integration_qag(&F, lower, upper, error_abs, error_rel, 100000,
                        //w, &result, &error); 
  gsl_integration_qag(&F, lower, upper, error_abs, error_rel, s,
                        2, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

double minimise(gsl_function &F, double lower, double upper, double guess, double error_abs, double error_rel){
  int status = 0, iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  double m = guess;
  double a = lower, b = upper;
  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc(T);
  gsl_min_fminimizer_set(s, &F, m, a, b);
  do{
    iter++;
    status = gsl_min_fminimizer_iterate(s);
    m = gsl_min_fminimizer_x_minimum(s);
    a = gsl_min_fminimizer_x_lower(s);
    b = gsl_min_fminimizer_x_upper(s);
    status = gsl_min_test_interval(a, b, 0.001, 0.0);
  } 
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_min_fminimizer_free(s);
  return m;


}
