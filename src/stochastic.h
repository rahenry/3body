#ifndef STOCHASTIC_H
#define STOCHASTIC_H

#include "general.h"
#include "three_body_system.h"

struct state_full;
struct characteristic_root_helper{
  characteristic_root_helper(three_body_system *tbs_, vector<double> *new_orthonormal_, vector<double> *h_) : tbs(tbs_), new_orthonormal(new_orthonormal_), h(h_) {};
  three_body_system *tbs;
  vector<double> *new_orthonormal;
  vector<double> *h;
};

struct random_state_test{
  state_full state;
  double eval;
  bool success;
  vector<double> S_column;
  vector<double> H_column;
};

double helper_function(double x, void * params);
bool competitive_sort_function(const random_state_test &t1, const random_state_test &t2);

#endif

