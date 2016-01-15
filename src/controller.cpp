#include "controller.h"
#include "utility.h"
#include "angular_subsystem.h"
#include <algorithm>
#include "three_body_system.h"

void controller::run_all(){
  // run all simulations specified by the various mass ratios, omega ratios, etc.

  n_simulations = (mass_ratios.size() > omega_ratios.size()) ? mass_ratios.size() : omega_ratios.size();
  n_simulations = (r0s_unlike.size() > n_simulations) ? r0s_unlike.size() : n_simulations;
  n_simulations = (r0s_like.size() > n_simulations) ? r0s_like.size() : n_simulations;

  print_vector(mass_ratios);

  for (unsigned int i = 0; i < n_simulations; i++){
    double omega_ratio = get_latest_element(omega_ratios, i);
    double mass_ratio = get_latest_element(mass_ratios, i);
    double r0_unlike = get_latest_element(r0s_unlike, i);
    double r0_like = get_latest_element(r0s_like, i);
    double delta_omega = (omega_ratio - 1.0) / (omega_ratio + 1.0);
    double delta_mass = (mass_ratio - 1.0) / (mass_ratio + 1.0);
    cout << delta_mass << " " << delta_omega << " " << r0_like << " " << r0_unlike << " " << endl;
    systems.push_back(
	three_body_system(this,
	  1.0 + 1.0 * delta_mass,
	  1.0 - delta_mass,
	  1.0 + 1.0 * delta_omega,
	  1.0 - delta_omega,
	  r0_like,
	  r0_unlike,
	  a0_like,
	  a0_unlike)
	);

    systems.back().competitive_loop();
    systems.back().calc_structure_factors_channels();

    if (calc_2body){
      cout << "Beginning two-body calculation..." << endl;
      systems.back().only_2body = 1;
      systems.back().basis_min_rho_scaled = basis_min_rho;
      systems.back().reset();
      systems.back().competitive_loop();

    }

  }
  


}

void controller::initialise(){
  mass_ratios = get_vector_double("mass_ratios", vector<double>(1, 1.0));
  omega_ratios = get_vector_double("omega_ratios", vector<double>(1, 1.0));

  permutation_number = get_int("permutation_number", 2);
  permutation_sign = get_int("permutation_sign", -1);
  cout << "permutaion_sign = " << permutation_sign << endl;
  parity = get_string("parity", "natural");
  r0s_unlike = get_vector_double("r0s_unlike", vector<double>(1, 0.003));
  r0s_like = get_vector_double("r0s_like", vector<double>(1, 0.003));
  a0_unlike = get_double("a0_unlike", 1.0E8);
  a0_like = get_double("a0_like", 1.0E8);
  use_interaction = get_bool("use_interaction", true);
  use_saved_states= get_bool("use_saved_states", true);

  ltotal = get_int("ltotal", 1);
  lmax_r = get_int("lmax_r", 3);
  lmax_rho = get_int("lmax_rho", 3);
  lmax_coup = get_int("lmax_coup", 3);
  lmax_com = get_int("lmax_com", 2);
  single_angular_state = get_bool("single_angular_state", 0);


  vector<int> max_states_default = vector<int>();
  max_states_default.push_back(100);
  max_states_default.push_back(150);
  max_states_default.push_back(150);
  max_states = get_vector_int("max_states", max_states_default);

  basis_max_r = get_double("basis_max_r", 2.00);
  basis_min_r = get_double("basis_min_r", 0.1);
  scale_min_to_r0_r = get_bool("scale_min_to_r0_r", true);

  basis_max_rho = get_double("basis_max_rho", 2.00);
  basis_min_rho = get_double("basis_min_rho", 0.1);
  scale_min_to_r0_rho = get_bool("scale_min_to_r0_rho", true);

  basis_max_com = get_double("basis_max_com", 2.00);
  basis_min_com = get_double("basis_min_com", 0.1);
  scale_min_to_r0_com = get_bool("scale_min_to_r0_com", false);

  n_eigen = get_int("n_eigen", 5);

  use_kinetic = get_bool("use_kinetic", true);
  use_trapping = get_bool("use_trapping", true);
  use_trapping_imbal = get_bool("use_trapping_imbal", use_trapping);
  calc_2body = get_bool("calc_2body", true);

  vector<int> sf_channels_default = vector<int>();
  sf_channels_default.push_back(2);
  sf_channels_default.push_back(3);
  sf_channels_default.push_back(0);
  sf_channels = get_vector_int("sf_channels", sf_channels_default);
  sf_minimum = get_double("sf_minimum", 0.0);
  sf_maximum = get_double("sf_maximum", 5.0);
  sf_n_points = get_int("sf_n_points", 200);

  use_spectrum = get_bool("use_spectrum", false);

  ainv_min = get_double("ainv_min", -10);
  ainv_max = get_double("ainv_max", 10);
  ainv_n_points = get_int("ainv_n_points", 10);

  state_inclusion_index = get_int("state_inclusion_index", -1);
  if (state_inclusion_index >= 0){
    lmax_r = state_inclusion_index;
    lmax_rho = state_inclusion_index;
    lmax_coup = state_inclusion_index;
    lmax_com = state_inclusion_index;
  }


  ang = new angular_subsystem(this);
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
