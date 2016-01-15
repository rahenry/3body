#ifndef CONTROLLER_H
#define CONTROLLER_H

struct controller{

  controller(string file_name) : name(file_name) {get_file_data(); initialise();};
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
  bool use_kinetic, use_trapping;

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

};

#endif
