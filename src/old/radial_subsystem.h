#ifndef RADIAL_SUBSYSTEM_H 
#define RADIAL_SUBSYSTEM_H 

#include "general.h"

struct three_body_system;
struct state_angular;
struct state_full;
struct angular_subsystem;

struct radial_subsystem{

  string name;

  radial_subsystem(three_body_system *tbs, angular_subsystem *as, double r0_);

  three_body_system *tbs;
  angular_subsystem *ang;

  // limits from angular subsystem
  int ltotal, lmax_r, lmax_rho, lmax_coup, lmax_com;
  bool single_angular_state, use_trapping_imbal;
  int n_states, n_states_E, n_states_transformed;
  int lmax, lim1, lim2, lim3, lim4, lim5, lim6, lim12, lim10p3, lim_kappa, lim_special, lim_cg, lim2_cg, m_offset_cg, lim_coup;
  string parity;

  // physical quantities
  double r0, mass_ratio, omega_ratio;
  double r0_like, r0_unlike, m1, m2, omega1, omega2;

  void generate_caches();
  void print_caches();

  vector<state_full> states;

  double scale_min_r, scale_max_r, scale_min_rho, scale_max_rho, scale_min_com, scale_max_com;

  uniform_real_distribution<double> rng_alpha_r, rng_alpha_rho, rng_alpha_com;
  state_full generate_new_state();
};


#endif
