#ifndef THREE_BODY_SYSTEM_H
#define THREE_BODY_SYSTEM

#include "states.h"
#include "caches.h"
#include "angular_subsystem.h"
#include <random>

struct three_body_system{

  three_body_system(string file_name) : name(file_name), masses{1.0, 1.0, 1.0}, omegas{1.0, 1.0, 1.0} {get_file_data(); initialise();};
  void initialise();

  // system name
  string name;

  // data from input file
  string file_data;

  // fundamental physical parameters
  vector<double> mass_ratios, omega_ratios;
  vector<double> r0s_unlike, r0s_like;
  double a0_unlike, a0_like;
  int permutation_number, permutation_sign;
  int ltotal;
  bool use_interaction;
  string parity;

  // fundamental numerical parameters
  double basis_max_r, basis_min_r, basis_max_rho, basis_min_rho, basis_max_com, basis_min_com;
  bool scale_min_to_r0_r, scale_min_to_r0_rho, scale_min_to_r0_com;
  int lmax_r, lmax_rho, lmax_com, lmax_coup;
  bool single_angular_state;

  // other options
  int n_eigen;
  bool use_kinetic, use_trapping, use_trapping_imbal, use_trapping_bal;

  // derived quantities and aliases
  double masses[3], omegas[3];
  double U[3][3], Uinv[3][3], trapping[3][3], kinetic[3][3];
  double m1, m2;
  double V0_unlike, V0_like, V0_prefactor;
  double r0_unlike, r0_like, r0_prefactor_unlike, r0_prefactor_like, r0m2_unlike, r0m2_like, r0_min;
  double *S, *H_full, *H_int, *H_non_int, *evals, *evecs;
  int n_runs;
  
  // larger caches
  vector<double> norms, basis_r, basis_rho, basis_com, gamma, integrals_gaussian_r, integrals_gaussian_r_int, integrals_gaussian_rho, integrals_gaussian_com;
  vector< vector< vector<double>* > > transformations;

  // functions to retrieve data from the input file
  void get_file_data();
  double get_double(string var_name, double default_value);
  int get_int(string var_name, int default_value);
  bool get_bool(string var_name, bool default_value);
  string get_string(string var_name, string default_value);
  vector<int> get_vector_int(string var_name, vector<int> default_name);
  vector<double> get_vector_double(string var_name, vector<double> default_value);

  // system control loops
  void loop_runs();
  // ***************************************
  // generation functions
  // ***************************************


  // coordinate matrices
  double U_element(int i, int j);
  double Uinv_element(int i, int j);
  double trapping_element(int i, int j);
  double kinetic_element(int i, int j);
  void generate_coordinate_matrices();
  void output_coordinate_matrices();
};

#endif
