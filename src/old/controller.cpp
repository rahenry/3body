#include "controller.h"
#include "utility.h"

void controller::initialise(){
  mass_ratios = get_vector_double("mass_ratios", vector<double>(1.0));
  omega_ratios = get_vector_double("omega_ratios", vector<double>(1.0));

  permutation_number = get_int("permutation_number", 1);
  permutation_sign = get_int("permutation_sign", -1);
  parity = get_string("parity", "all");
  r0s_unlike = get_vector_double("r0s_unlike", vector<double>(0.01));
  r0s_like = get_vector_double("r0s_like", vector<double>(0.01));
  a0_unlike = get_double("a0_unlike", 1.0E8);
  a0_like = get_double("a0_like", 1.0E8);
  use_interaction = get_bool("use_interaction", true);

  ltotal = get_int("ltotal", 1);
  lmax_r = get_int("lmax_r", 2);
  lmax_rho = get_int("lmax_rho", 1);
  lmax_coup = get_int("lmax_coup", 1);
  lmax_com = get_int("lmax_com", 1);
  single_angular_state = get_bool("single_angular_state", 0);

  basis_max_r = get_double("basis_max_r", 1.10);
  basis_min_r = get_double("basis_min_r", 0.8);
  scale_min_to_r0_r = get_bool("scale_min_to_r0_r", true);

  basis_max_rho = get_double("basis_max_rho", 1.10);
  basis_min_rho = get_double("basis_min_rho", 0.8);
  scale_min_to_r0_rho = get_bool("scale_min_to_r0_rho", true);

  basis_max_com = get_double("basis_max_com", 1.10);
  basis_min_com = get_double("basis_min_com", 0.8);
  scale_min_to_r0_com = get_bool("scale_min_to_r0_com", false);

  n_eigen = get_int("n_eigen", 5);

  use_kinetic = get_bool("use_kinetic", true);
  use_trapping = get_bool("use_trapping", true);
}


void controller::get_file_data(){
  ifstream input_stream(name);
  stringstream buffer;
  buffer << input_stream.rdbuf();
  file_data = buffer.str();
}

int controller::get_int(string var_name, int default_value){
  return extract_int(file_data, var_name, default_value);
}

double controller::get_double(string var_name, double default_value){
  return extract_double(file_data, var_name, default_value);
}

bool controller::get_bool(string var_name, bool default_value){
  return extract_bool(file_data, var_name, default_value);
}

string controller::get_string(string var_name, string default_value){
  return extract_string(file_data, var_name, default_value);
}

vector<double> controller::get_vector_double(string var_name, vector<double> default_value){
  return extract_vector_double(file_data, var_name, default_value);
}

vector<int> controller::get_vector_int(string var_name, vector<int> default_value){
  return extract_vector_int(file_data, var_name, default_value);
}
