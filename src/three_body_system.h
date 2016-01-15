#ifndef THREE_BODY_SYSTEM_H
#define THREE_BODY_SYSTEM_H

#include "general.h"
#include "states.h"
#include "gsl/gsl_blas.h"
#include "B_coef_cache.h"

struct controller;
struct angular_subsystem;
struct random_state_test;

struct three_body_system{

  three_body_system(controller *cont_, double m1_, double m2_, double omega1_, double omega2_, double r0_like_, double r0_unlike_, double a0_like_, double a0_unlike_);
  void initialise();

  // system name
  string name;

  // pointers to other components
  controller *cont;
  angular_subsystem *ang;

  bool bb;
  int cc;

  // fundamental physical parameters
  double m1, m2;
  double omega1, omega2;
  double r0_like, r0_unlike;
  double a0_like, a0_unlike;
  bool only_2body;

  int lim1, lim_tr, lim_kappa;

  vector<int> channels_allowed;

  // derived quantities and aliases
  double masses[3], omegas[3];
  double trapping[6][3][3], kinetic[6][3][3];
  gsl_matrix *U_matrix_base, *Uinv_matrix_base;
  double scale_factor_r[6], scale_factor_rho[6], scale_factor_com[6];
  vector< gsl_matrix* > U, Uinv;
  double V0_unlike, V0_like;
  double r0_prefactor_unlike, r0_prefactor_like, r0m2_unlike, r0m2_like, r0_min;
  vector<double> S, H, evals, evecs;
  int n_runs;
  int basis_size;
  double basis_max_r, basis_max_rho, basis_max_com, basis_min_r_scaled, basis_min_rho_scaled, basis_min_com_scaled;
  
  vector<double> gamma, gamma2, fact, fact2;
  void generate_factorials();
  vector< vector<double> > single_particle_transformations;
  vector<state_full> states;

  double U_element(int i, int j);
  double Uinv_element(int i, int j);
  double trapping_element(int i, int j, int channel);
  double kinetic_element(int i, int j, int channel);
  void generate_coordinate_matrices();
  void output_coordinate_matrices();
  int channel_index(int channel);
  void report_length_scales();


  void initialise_transformations();
  void get_single_particle_transformation(vector<double> &tar, int channel);
  void calc_transformation(vector<double> &tar, state_full &s);
  void calc_transformation_coefficients(vector<double> &tar, vector<double> &t);
  void calc_transformed_sums(vector<double> &tar, state_full &s1, state_full &s2);
  void apply_transformation(vector<double> &tar, state_full &s1, state_full &s2, vector<double> &tr1, vector<double> &tr2);

  void get_permutation_matrix(vector<double> &tar, int channel);
  void calc_coordinate_transformation(vector<double> &tar, int channel_initial, int channel_final);
  void calc_basis_transformation(vector<double> &tar, vector<double> &t);
  void calc_transformed_state(vector<double> &tar, state_full &s1, state_full &s2, int channel_final);

  vector< vector<B_coef_cache> > B_coef_caches;
  void generate_B_coef_caches();
  double B(state_full &state_start, int transformation_index, int channel_end);

  int get_kappa_limit(state_angular *a1, state_angular *a2);
  double gaussian_integral(double a, int b);
  double bessel_integral(vector<double> &transformed_state, int kappa, int nr, int nrho, double interaction);
  void calc_matrix_elements(vector<double> &tar, state_full s1, state_full s2);
  void loop_transformations_H(vector<double> &tar, state_full &s1, state_full &s2, int channel_working, double prefactor, bool calc_interaction_only);
  void calc_matrix_element_contributions(vector<double> &tar, state_full &s1, state_full &s2, int transformation_index1, int transformation_index2, double prefactor, int channel_working, bool calc_interaction_only);
  void calc_angular_couplings(vector<double> &tar, state_angular *a1, state_angular *a2);

  // state generation and manipulation
  void generate_states_geometric();
  void geometric_series(vector<double> &tar, double geo_min, double geo_max, double N);
  void reset();

  // solution methods
  void regenerate_matrices();
  void regenerate_norms();
  void diagonalise();


  // stochastic methods
  bool focus_states_near_r0;
  void stochastic_loop();
  state_full generate_random_state(int state_generation_mode);
  state_full generate_random_state(double min_r, double max_r, double min_rho, double max_rho, double min_com, double max_com, int angular_index, int channel);
  state_full generate_refining_state(state_full &s);
  void calc_new_matrix_elements(state_full &new_state, vector<double> &S_new_column, vector<double> &H_new_column);
  void calc_eigenstate_matrix_elements(state_full &new_state, vector<double> &s, vector<double> &S_new_column, vector<double> &H_new_column);
  void calc_orthonormal_test_state(vector<double> &tar, state_full &new_state, vector<double> &S_new_column, vector<double> &s);
  void calc_orthonormal_hamiltonian_elements(vector<double> &h, vector<double> &new_orthonormal, vector<double> &S_new_column, vector<double> &H_new_column);
  double characteristic_function(double E, vector<double> &new_orthonormal, vector<double> &h);
  bool find_characteristic_root(double &tar, vector<double> &new_orthonormal, vector<double> &h, int N);
  void test_state(random_state_test &t);
  void adjoin_new_columns(vector<double> &S_new_column, vector<double> &H_new_column);

  // competitive method
  //
  bool select_new_state_competitive(int n_states, int state_generation_mode);
  void competitive_loop();
  void competitive_state_loop();
  bool competitive_end_test(double eval_old, double eval_new);


  //tests
  void evecs_test();
  void new_orthonormal_test(vector<double> &new_orthonormal, vector<double> &S_new_column);
  double previous_eigenvalue_change;

  // structure factors
  vector< vector<double> > sf_results;
  vector<double> sf_radii;
  double calc_sf_channel_loop_transformations(double rP, state_full &s1, state_full &s2, int channel_working, double prefactor);
  void calc_structure_factors_channels();
  double calc_sf_channel_element(double rP, state_full s1, state_full s2, int channel_working);
  double calc_sf_channel_contribution(double rP, state_full &s1, state_full &s2, int transformation_index1, int transformation_index2, int channel_working, double prefactor);
  double calc_sf_channel_summand(int kappa, int nr, int nrho, vector<double> transformed_state, double rP);
  void print_structure_factors();

  // hyperspherical
  double hyper_rescale_param;
  double length_hyperspherical;
  double calc_sf_hyper_summand(int kappa, int nr, int lr, int nrho, int lrho, vector<double> transformed_state, double h);

  // testing of structure factors
  bool test_structure_factors;
  vector< vector<double> > sf_test_matrices;
  vector<double> sf_test_results;
  void compute_sf_test_results();

  // spectrum
  void generate_spectrum();

  friend ostream &operator<<(ostream &os, three_body_system const &t);
};


double epsilon(int kappa, double w);
#endif
