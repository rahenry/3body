#include "three_body_system.h"
#include "angular_subsystem.h"
#include "states.h"
#include "gsl/gsl_blas.h"
#include "controller.h"
#include "utility.h"
#include "eigenvalues.h"
#include "rng.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "stochastic.h"
#include <tuple>
#include <algorithm>
#include <omp.h>

void three_body_system::new_orthonormal_test(vector<double> &new_orthonormal, vector<double> &S_new_column){
  bool flag = false;
  for (int i = 0; i < basis_size; i++){
    double overlap = 0.0;
    for (int a = 0; a < basis_size; a++){
      for (int b = 0; b < basis_size+1; b++){
	double x1, x2;
	x1 = evecs[a + i * basis_size];
	x2 = new_orthonormal[b];
	double S_elem;
	if (b == basis_size) S_elem = S_new_column[a];
	else S_elem = S[a+basis_size*b];
	overlap += x1 * x2 * S_elem;
      }
    }
    cout << "Overlap < " << i << " | new " << " > = " << overlap << endl;
  }
}

void three_body_system::evecs_test(){
  bool flag = false;
  for (int i = 0; i < basis_size; i++){
    for (int j = 0; j < basis_size; j++){
      double overlap = 0.0;
      for (int a = 0; a < basis_size; a++){
	for (int b = 0; b < basis_size; b++){
	  double x1, x2;
	  if (flag){
	    x1 = evecs[i + a * basis_size];
	    x2 = evecs[j + b * basis_size];
	  }
	  else{
	    x1 = evecs[a + i * basis_size];
	    x2 = evecs[b + j * basis_size];
	  }


	  overlap += x1 * x2 * S[a + basis_size * b];
	}
      }
      cout << "Overlap < " << i << " | " << j << " > = " << overlap << endl;
    }
  }


}


void three_body_system::competitive_loop(){

  int channel_0 = 3;
  int cind = channel_index(channel_0);

  if (cont->use_saved_states){
    string state_input_name = name + ".states";
    if (only_2body) state_input_name += "_2body";
    ifstream state_input(state_input_name);
    string line;

    while (getline(state_input, line)){
      string line_temp = line;
      if (line_temp.find("Fin") == std::string::npos)
        states.push_back(state_full(line, ang));
    }

    cout << "Successfully read states from file." << endl;

  }
  if (states.size() == 0){
    states.push_back(
      state_full(0, &ang->states[0], 
	scale_factor_r[cind], scale_factor_rho[cind], scale_factor_com[cind],
	channel_0)
      );
    string state_output_name = name + ".states";
    if (only_2body) state_output_name += "_2body";
    ofstream state_output(state_output_name, ios_base::app);
    state_output << states.back();
  }

  regenerate_matrices();
  diagonalise();


  // determine what mode we're in initially...
  int current_mode = 0;
  int state_counter = states.size();
  for (auto mode_size : cont->max_states){
    if (state_counter >= mode_size){ 
      current_mode++;
      state_counter -= mode_size;
    }
    else{
      break;
    }

  }

  int end_test_counter = 0;
  int fail_counter = 0;
  int fail_limit = 2;

  double energy_old = evals[0];
  double convtol = 1E-6;

  string energy_output_name = name + ".energies";
  if (only_2body){
    energy_output_name = name + ".2body";
  }
  ofstream energy_output(energy_output_name, ios_base::app);
  for (; current_mode < int(cont->max_states.size()); current_mode++){
    int mode_size = cont->max_states[current_mode];

    int fail_counter = 0;
    cout << "sc = " << state_counter << endl;
    while (state_counter++ <= mode_size){
      if (select_new_state_competitive(1000, current_mode)){
        if (relative_difference(energy_old, evals[0]) < convtol){
          fail_counter++;
        }
        else{
          cout << energy_old << " " << evals[0] << " " << relative_difference(energy_old, evals[0]) << endl;
          energy_old = evals[0];
          fail_counter = 0;
          energy_output << states.size() << ": " << flush;

          for (unsigned int i=0; i<10; i++){
            if (i < evals.size()){
              energy_output << evals[i] << " ";
            }
          }
          energy_output << endl;
          //break;
        }
      }
      else fail_counter++;
      if (fail_counter != 0) cout << "fails = " << fail_counter << endl;
      if (fail_counter >= fail_limit){
        cout << "Cannot find new states, ending simulation..." << endl;
        break;
      }
    }
    state_counter = 0;
    cout << "Completed mode " << current_mode << endl;
    cout << "state_counter = " << state_counter << endl;
    cout << "mode_size = " << mode_size << endl;
  }

  //print_matrix_to_file(&S[0], basis_size, name + ".S");
  //print_matrix_to_file(&H[0], basis_size, name + ".H");

}

bool three_body_system::competitive_end_test(double eval_old, double eval_new){
  double diff_abs = abs(eval_old - eval_new);
  double diff_rel = relative_difference(eval_old, eval_new);

  return (diff_abs < 0.0001 && diff_rel < 0.0001);
}

bool competitive_sort_function(const random_state_test &t1, const random_state_test &t2){

  if (t1.success && !t2.success) return true;
  if (t2.success && !t1.success) return false;
  else return (t1.eval < t2.eval);
}

bool three_body_system::select_new_state_competitive(int n_states, int state_generation_mode){
  int quiet = 0;

  clock_t search_time_start = clock();

  bool found_state_success = false;
  
  vector<random_state_test> results;
  results.resize(n_states);

  omp_set_num_threads(omp_get_num_threads() - 1);
  #pragma omp parallel for
  for (int k=0; k < n_states; k++){
    results[k] = random_state_test();
    results[k].state = generate_random_state(state_generation_mode);
    test_state(results[k]);
  }

  sort(results.begin(), results.end(), competitive_sort_function);

  for (int k=0; k < n_states; k++){
    if (results[k].success) found_state_success = true;
  }

  if (!found_state_success){
    cout << "Search failed, t = " << double(clock() - search_time_start) / (double)CLOCKS_PER_SEC << endl;
    return false;
  }
  if (found_state_success){
    int refining_counter = 0;
    vector<random_state_test> results_refining;
    while(refining_counter++ < 50){
      results_refining.push_back(random_state_test());
      results_refining.back().state = generate_refining_state(results[0].state);
      //cout << "generated refining state: " << results_refining.back().state;
      test_state(results_refining.back());
    }

    sort(results_refining.begin(), results_refining.end(), competitive_sort_function);
    bool found_refining_state_success = false;
    //cout << "refining info: " << results_refining.size() << endl;
    for (auto &r : results_refining){
      if (r.success) found_refining_state_success = true;
    }
    //cout << "found_refining_state_success = " << found_refining_state_success << endl;

    if (found_refining_state_success){
      states.push_back(results_refining[0].state);
      states.back().index = basis_size++;
      adjoin_new_columns(results_refining[0].S_column, results_refining[0].H_column);
    }
    else{
      states.push_back(results[0].state);
      states.back().index = basis_size++;
      adjoin_new_columns(results[0].S_column, results[0].H_column);
    }


    double search_time = double(clock() - search_time_start) / (double)CLOCKS_PER_SEC;

    clock_t diag_time_start = clock();
    diagonalise();
    double diag_time = double(clock() - diag_time_start) / (double)CLOCKS_PER_SEC;

    cout << "t = " << search_time << " " << diag_time << endl;


    string state_output_name = name + ".states";
    if (only_2body) state_output_name += "_2body";
    ofstream state_output(state_output_name, ios_base::app);
    state_output << states.back();

  }

  if (found_state_success && !quiet) cout << states.back();
  return found_state_success;

}

void three_body_system::adjoin_new_columns(vector<double> &S_new_column, vector<double> &H_new_column){
  vector<double> S_new, H_new;

  for (int i=0; i<basis_size; i++){
    for (int j=0; j<basis_size; j++){
      if (i == basis_size - 1){
	S_new.push_back(S_new_column[j]);
	H_new.push_back(H_new_column[j]);
      }
      else if (j == basis_size - 1){
	S_new.push_back(S_new_column[i]);
	H_new.push_back(H_new_column[i]);
      }
      else{
	int index = i + (basis_size - 1) * j;
	S_new.push_back(S[index]);
	H_new.push_back(H[index]);
      }
    }
  }

  S = S_new;
  H = H_new;
}

void three_body_system::test_state(random_state_test &t){

  calc_new_matrix_elements(t.state, t.S_column, t.H_column);
  vector<double> s, h;
  calc_eigenstate_matrix_elements(t.state, s, t.S_column, t.H_column);
  vector<double> new_orthonormal;
  calc_orthonormal_test_state(new_orthonormal, t.state, t.S_column, s);
  //new_orthonormal_test(new_orthonormal, t.S_column);
  calc_orthonormal_hamiltonian_elements(h, new_orthonormal, t.S_column, t.H_column);
  t.success = find_characteristic_root(t.eval, new_orthonormal, h, 0);
}

void three_body_system::stochastic_loop(){

  // start with a sigle state, the QHO ground state
  int channel_0 = 3;
  int cind = channel_index(channel_0);
  states.push_back(
      state_full(0, &ang->states[0], 
	scale_factor_r[cind], scale_factor_rho[cind], scale_factor_com[cind],
	channel_0)
      );
  if (true){
    states.push_back(
	state_full(1, &ang->states[0], 
	  0.8*scale_factor_r[cind], 0.8*scale_factor_rho[cind], scale_factor_com[cind],
	  channel_0)
	);
  }
  // create initial matrix and "diagonalise"

  regenerate_matrices();
  diagonalise();
  select_new_state_competitive(10, 0);


    /*state_full new_state = generate_random_state();
    cout << "New state... " << new_state;

    double eval;
    vector<double> S_new_column, H_new_column;

    test_state(new_state, S_new_column, H_new_column, eval);



    cout << "New eigenvalue = " << eval << endl;*/


}

double helper_function(double x, void * params){
  characteristic_root_helper *helper = (characteristic_root_helper *) params;
  return helper->tbs->characteristic_function(x, *helper->new_orthonormal, *helper->h);
}

bool three_body_system::find_characteristic_root(double &tar, vector<double> &new_orthonormal, vector<double> &h, int N){
  // find the Nth root of the characteristic function specified by new_orthonormal and h

  gsl_set_error_handler_off();

  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  gsl_function F;

  F.function = &helper_function;
  characteristic_root_helper * params = new characteristic_root_helper(this, &new_orthonormal, &h);
  F.params = params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);

  double x_lo = 0.0;
  double x_hi = evals[N] - 0.00001;
  if (N > 0) x_lo = evals[N-1]*1.000000001;
  else{
    x_lo = x_hi - 1.0;
    while (true){
      double temp = GSL_FN_EVAL(&F, x_lo);
      //cout << x_lo << " " << temp << endl;
      if (temp < 0.0) break;
      x_lo -= 1.0;
      if (x_lo < x_hi - 10.0) break;

    }
    //if (evals[0] < 0.0) x_lo = evals[0] * 1.5;
    //else x_lo = min(0.0, evals[0] - 250.0);
  }


  status = gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    //cout << x_lo << " " << x_hi << endl;
    //cout << GSL_FN_EVAL(&F, x_lo) << " " << GSL_FN_EVAL(&F, x_hi) << endl;
  if (status != GSL_SUCCESS){
    //cout << "Root solver failure?" << endl;
    return false;
  }

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,
                                       0, 1E-8);

    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  tar = r;

  return true;
}


void three_body_system::calc_orthonormal_hamiltonian_elements(vector<double> &h, vector<double> &new_orthonormal, vector<double> &S_new_column, vector<double> &H_new_column){
  h = vector<double>();
  h.reserve(basis_size);


  for (int k=0; k<basis_size; k++){
    h.push_back(0.0);
    for (int i=0; i<basis_size; i++){
      int index = i + basis_size * k;
      for (int j=0; j<basis_size+1; j++){
	double h_elem;
	if (j < basis_size) h_elem = H[i + j*basis_size];
	else h_elem = H_new_column[i];
	h.back() += evecs[index] * new_orthonormal[j] * h_elem;
      }
    }
  }

  h.push_back(0.0);
  for (int i=0; i<basis_size+1; i++){
    for (int j=0; j<basis_size+1; j++){
      double h_elem;
      if (i < basis_size && j < basis_size) h_elem = H[i + j*basis_size];
      else if (i < basis_size) h_elem = H_new_column[i];
      else if (j < basis_size) h_elem = H_new_column[j];
      else h_elem = H_new_column.back();

      double temp =  new_orthonormal[i] * new_orthonormal[j] * h_elem;
      h.back() += temp;
    }
  }

}

double three_body_system::characteristic_function(double E, vector<double> &new_orthonormal, vector<double> &h){
  double sum1 = 0.0;
  for (int i=0; i<basis_size; i++){
    sum1 += pow(h[i], 2) / (E - evals[i]);
  }

  double temp = E - h.back() - sum1;
  return temp;


}

void three_body_system::calc_orthonormal_test_state(vector<double> &tar, state_full &new_state, vector<double> &S_new_column, vector<double> &s){
  double sum1 = 0.0;
  for (int i=0; i<basis_size; i++){
    sum1 += pow(s[i], 2);
    //double overlap = 
      //gaussian_integral(states[i].alpha_r+new_state.alpha_r, 0) *
      //gaussian_integral(states[i].alpha_rho+new_state.alpha_rho, 0) *
      //gaussian_integral(states[i].alpha_com+new_state.alpha_com, 0);
    //sum1 += pow(states[i].norm, 2);
    //sum1 += overlap*overlap;
  }
  double norm = pow(
      (S_new_column.back() - sum1),
      -0.5
      );

  tar = vector<double>(basis_size+1, 0.0);
  tar.back() = norm * pow(new_state.norm, -2);
  tar.back() = norm * pow(S_new_column.back(), -0.5);
  tar.back() = norm * pow(new_state.norm, -2);
  tar.back() = norm;

  

  for (int i=0; i<basis_size; i++){
    for (int j=0; j<basis_size; j++){
      int index = j + basis_size * i;
      tar[j] -= norm * evecs[index] * s[i];
    }
  }
}

void three_body_system::calc_eigenstate_matrix_elements(state_full &new_state, vector<double> &s, vector<double> &S_new_column, vector<double> &H_new_column){
  s = vector<double>();
  s.reserve(basis_size);

  for (int i=0; i<basis_size; i++){
    s.push_back(0.0);
    for (int j=0; j<basis_size; j++){
      int index = j + i * basis_size;
      s.back() += evecs[index] * S_new_column[j];
    }
  }
}

void three_body_system::calc_new_matrix_elements(state_full &new_state, vector<double> &S_new_column, vector<double> &H_new_column){

  S_new_column.reserve(basis_size);
  H_new_column.reserve(basis_size);

  vector<double> contribs;
  double H_elem;
  for (auto & s : states){
    calc_matrix_elements(contribs, s, new_state);
    S_new_column.push_back(contribs[0]);
    H_elem = 0.0;
    for (unsigned int i=1; i<contribs.size(); i++) H_elem += contribs[i];
    H_new_column.push_back(H_elem);
  }
  calc_matrix_elements(contribs, new_state, new_state);
  S_new_column.push_back(contribs[0]);
  H_elem = 0.0;
  for (unsigned int i=1; i<contribs.size(); i++) H_elem += contribs[i];
  H_new_column.push_back(H_elem);

}

state_full three_body_system::generate_random_state(int state_generation_mode){
  int quiet = 1;
  int channel = 0;
  while (abs(channel) == 0){
    channel = random_int(-3, 3);
  }
  channel = 3;
  int cind = channel_index(channel);

  double random_log_r = random_double( log(basis_min_r_scaled), log(basis_max_r) );
  if (true){
    int randomiser = random_int(0,3);
    if (randomiser != 0){
      random_log_r = random_double( log(0.9*basis_min_r_scaled), log(1000.0 * basis_min_r_scaled) );
    }


  }
  double alpha_r = scale_factor_r[cind] * exp(random_log_r);
  double random_log_rho = random_double( log(basis_min_rho_scaled), log(basis_max_rho) );
  double alpha_rho = scale_factor_rho[cind] * exp(random_log_rho);
  double alpha_com = scale_factor_com[cind] * random_double(basis_min_com_scaled, basis_max_com);
  
  int temp_int = -1 + ang->generation_modes.size();
  int temp_mode = min(state_generation_mode, temp_int);
  int temp_mode_max = ang->generation_modes[temp_mode].size() - 1;
  int angular_index = ang->generation_modes[temp_mode][random_int(0, temp_mode_max)];
  if (only_2body){
    while (true){
      angular_index = ang->generation_modes[temp_mode][random_int(0, temp_mode_max)];
      if (&ang->states[angular_index].lcom == 0
	  && ang->states[angular_index].lcoup == ang->ltotal){
	break;
      }
    }
  }

  if (!quiet){
    cout << "Generated random state: " << state_full(-1, &ang->states[angular_index], alpha_r, alpha_rho, alpha_com, channel);
  }

  return state_full(states.size(), &ang->states[angular_index], pow(alpha_r,-2), pow(alpha_rho,-2), pow(alpha_com,-2), channel);
}

state_full three_body_system::generate_random_state(double min_r, double max_r, double min_rho, double max_rho, double min_com, double max_com, int angular_index, int channel){
  int quiet = 1;
  int cind = channel_index(channel);

  double random_log_r = random_double( log(min_r), log(max_r) );
  double alpha_r = scale_factor_r[cind] * exp(random_log_r);
  double random_log_rho = random_double( log(min_rho), log(max_rho) );
  double alpha_rho = scale_factor_rho[cind] * exp(random_log_rho);
  double alpha_com = scale_factor_com[cind] * random_double(min_com, max_com);
  
  return state_full(states.size(), &ang->states[angular_index], pow(alpha_r,-2), pow(alpha_rho,-2), pow(alpha_com,-2), channel);
}

state_full three_body_system::generate_refining_state(state_full &s){
  return state_full(s.index,
      s.ang,
      random_double(0.8*s.alpha_r, 1.2*s.alpha_r),
      random_double(0.9*s.alpha_rho, 1.1*s.alpha_rho),
      random_double(1.0*s.alpha_com, 1.0*s.alpha_com),
      s.channel);
}

