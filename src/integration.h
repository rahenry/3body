#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "general.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>

double integrate(gsl_function &F, double lower, double upper, double error_abs, double error_rel);
double integrate_to_inf(gsl_function &F, double lower, double error_abs, double error_rel);
double minimise(gsl_function &F, double lower, double upper, double guess, double error_abs, double error_rel);


#endif
