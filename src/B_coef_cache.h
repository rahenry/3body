#ifndef B_COEF_CACHE_H
#define B_COEF_CACHE_H

#include "general.h"
#include "states.h"

struct three_body_system;
struct B_coef_cache{

  B_coef_cache(three_body_system *tbs, int channel_start_, int channel_end_);

  three_body_system* tbs;

  int channel_start, channel_end;
  vector< vector<double> > B_coefs;

  double calc_B_coef(state_full &state_start, state_full &state_end);

  double T11, T12, T21, T22;
  double D(int L, int k1, int l1, int k2, int l2);


};

#endif
