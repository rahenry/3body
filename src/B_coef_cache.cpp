#include "B_coef_cache.h"
#include "angular_subsystem.h"
#include "three_body_system.h"
#include "states.h"
#include "utility.h"

B_coef_cache::B_coef_cache(three_body_system *tbs_, int channel_start_, int channel_end_) : tbs(tbs_), channel_start(channel_start_), channel_end(channel_end_){

  bool quiet = 1;

  // we permute the start state and change end_channel to 3, which is equivalent

  B_coefs = vector< vector<double> >();

  state_full dummy_state_start(0, &tbs->ang->states[0], 1.0, 1.0, 1.0, channel_start);
  state_full dummy_state_end(0, &tbs->ang->states[0], 1.0, 1.0, 1.0, channel_end);

  vector<double> transformation;
  tbs->calc_coordinate_transformation(transformation, channel_start, channel_end);
  tbs->calc_coordinate_transformation(transformation, channel_end, channel_start);
  T11 = transformation[0];
  T12 = transformation[1];
  T21 = transformation[3];
  T22 = transformation[4];

  if (!quiet) cout << "Generating a B cache, " << channel_start << setw(2) << " -> " << channel_end << setw(2) << ", T = {" << T11 << ", " << T12 << "; " << T21 << ", " << T22 << "}" << endl;

  for (auto sa : tbs->ang->states){
    B_coefs.push_back(vector<double>());
    state_full state_start = state_full(dummy_state_start, &sa);


    for (auto transformation_index : state_start.ang->allowed_transformations){
      state_full state_end = state_full(dummy_state_end, &tbs->ang->states_transformed[transformation_index]);

      B_coefs.back().push_back(calc_B_coef(state_start, state_end));


    }

  }


}

double B_coef_cache::calc_B_coef(state_full &state_start, state_full &state_end){
  bool quiet = 1;
  if (channel_end == -1) quiet = true;
  if (channel_start != 3) quiet = true;
  state_angular *s = state_start.ang;
  state_angular *e = state_end.ang;

  int lim1 = s->lr + 2*s->nr;
  int lim2 = s->lrho + 2*s->nrho;
  int lim3 = e->lr + 2*e->nr;
  int lim4 = e->lrho + 2*e->nrho;
  int temp3 = 0;
  int lim_r1 = min(lim1, lim3) + temp3;
  int lim_r2 = min(lim1, lim4) + temp3;
  int lim_rho1 = min(lim2, lim3) + temp3;
  int lim_rho2 = min(lim2, lim4) + temp3;

  if (!quiet) cout << state_start.ang->index << " " << state_end.ang->index << ": ";

  double res = 0.0;
  for (int lr1=0; lr1 <= lim_r1; lr1++)
  for (int nr1=0; nr1 <= lim_r1/2; nr1++)
  for (int lrho1=0; lrho1 <= lim_rho1; lrho1++)
  for (int nrho1=0; nrho1 <= lim_rho1/2; nrho1++)
  for (int lr2=0; lr2 <= lim_r2; lr2++)
  for (int nr2=0; nr2 <= lim_r2/2; nr2++)
  for (int lrho2=0; lrho2 <= lim_rho2; lrho2++)
  for (int nrho2=0; nrho2 <= lim_rho2/2; nrho2++){
    if (2*nr1 + lr1 + 2*nr2 + lr2 != lim1 ||
        2*nrho1 + lrho1 + 2*nrho2 + lrho2 != lim2 ||
        2*nr1 + lr1 + 2*nrho1 + lrho1 != lim3 ||
        2*nr2 + lr2 + 2*nrho2 + lrho2 != lim4)
      continue;

    if (!triangle(lr1, lr2, s->lr) ||
        !triangle(lrho1, lrho2, s->lrho) ||
        !triangle(lr1, lrho1, e->lr) ||
        !triangle(lr2, lrho2, e->lrho))
      continue;

    double E = tbs->ang->E(lr1, lr2, s->lr, lrho1, lrho2, s->lrho, e->lr, e->lrho, s->lcoup);
    if (E == 0.0) continue;

    double temp1 = (T12 == 0.0 && 2*nr2 + lr2 == 0) ? 1.0 : pow(T12, 2*nr2 + lr2);
    double temp2 = (T21 == 0.0 && 2*nrho1 + lrho1 == 0) ? 1.0 : pow(T21, 2*nrho1 + lrho1);
    double D1 = D(s->lr, nr1, lr1, nr2, lr2);
    double D2 = D(s->lrho, nrho1, lrho1, nrho2, lrho2);
    res += E * D1 * D2 * pow(T11, 2*nr1 + lr1) * pow(T22, 2*nrho2 + lrho2) * temp1 * temp2;
    if (!quiet && (state_start.ang->index == 1 && state_end.ang->index == 2)) cout << "! " << E << " " << temp1 << " " << temp2 << " " << D1 << " " << D2 << " | " << lr1 << " " << lr2 << " " << lrho1 << " " << lrho2 << ", lims = " << lim1 << " " << lim2 << " " << lim3 << " " << lim4 << endl;

  }

  if (!quiet) cout << "result: " << res << endl;


  return res;


}

double B_coef_cache::D(int L, int k1, int l1, int k2, int l2){
  return 2.0 * sqrt(Pi*(2.0*l1+1.0)*(2.0*l2+1.0)/(2.0*L+1.0))
          * tbs->ang->cg(l1,0,l2,0,L,0) * tbs->fact2[2*L+1]
          / (pow(2.0,k1+k2) * tbs->fact[k1] * tbs->fact[k2] * tbs->fact2[2*k1+2*l1+1]
          * tbs->fact2[2*k2+2*l2+1]);

}
