#ifndef ANGULAR_SUBSYSTEM_H
#define ANGULAR_SUBSYSTEM_H

#include "general.h"

struct controller;
struct state_angular;
struct angular_subsystem{

  string name;

  angular_subsystem(controller *cont_);

  controller *cont;

  int ltotal, lmax_r, lmax_rho, lmax_coup, lmax_com;
  bool single_angular_state, use_trapping_imbal;
  int n_states, n_states_E, n_states_transformed;
  int lmax, lim1, lim2, lim3, lim4, lim5, lim6, lim12, lim10p3, lim_tr, lim_kappa, lim_special, lim_cg, lim2_cg, m_offset_cg, lim_coup;
  string parity;

  vector<state_angular> states, states_transformed;
  void generate_states();
  void generate_states_transformed();
  void print_angular_states();

  vector<double> cache_cg, cache_E, cache_D, cache_P, cache_Q, cache_F, cache_G_r, cache_G_rho, cache_H_r, cache_H_rho;

  double cg(int l1, int m1, int l2, int m2, int l3, int m3);
  double E(int l1, int l2, int l3, int l4, int l5, int l6, int l7, int l8, int l9);
  double P(int l1, int m1, int l2, int m2, int l3, int m3, int l4, int m4);
  double Q(state_angular *s1, state_angular *s2, int a, int b, int mcom);
  double F(state_angular *s1, state_angular *s2, int kappa);
  double G_r(state_angular *s1, state_angular *s2, int a);
  double G_rho(state_angular *s1, state_angular *s2, int b);
  double H_r(state_angular *s1, state_angular *s2, int kappa);
  double H_rho(state_angular *s1, state_angular *s2, int kappa);

  vector<vector<int> > generation_modes;
  void generate_generation_modes();
  bool generation_mode_test(state_angular &s, int mode);

  void generate_cg();
  void generate_E();
  void generate_P();
  void generate_Q();
  void generate_F();
  void generate_G();
  void generate_H();


  void generate_caches();
  void print_caches();

  uniform_int_distribution<int> rng_angular;
  state_angular* random_state();
};


#endif
