#include "general.h"
#include "tests.h"
#include "controller.h"
#include "three_body_system.h"
#include "angular_subsystem.h"
#include "utility.h"
#include "states.h"
#include "rng.h"
#include "angular_subsystem.h"
#include "stochastic.h"
#include <gsl/gsl_sf_bessel.h>
#include "integration.h"
#include "structure_factors.h"

void test6(){
  controller cont1("test6");
  three_body_system tbs(&cont1, 1.0, 1.0, 2.0, 1.0, 0.01, 0.01, 1E-8, 1E-8);
  cont1.ang->print_angular_states();
  tbs.focus_states_near_r0 = true;
  tbs.generate_states_geometric();
  tbs.regenerate_matrices();
  tbs.diagonalise();
  int i1 = 1; int i2 = 0;
  cout << tbs.states[i1] << endl << tbs.states[i2] << endl;
  vector<double> res1;
  tbs.calc_matrix_elements(res1, tbs.states[i1], tbs.states[i2]);
  print_vector(res1);
  cont1.ang->print_caches();

  //testP(cont1.ang);
}

void test5(){
  controller cont1("test5");
  cout << "permutation_number = " << cont1.permutation_number << endl;
  cout << "permutation_sign = " << cont1.permutation_sign << endl;
  three_body_system tbs(&cont1, 12.313, 1.00, 1.0, 3.0, 0.01, 0.01, 1E-8, 1E-8);
  tbs.focus_states_near_r0 = true;
  tbs.competitive_loop();
  //tbs.stochastic_loop();
}


void test4(){
  cout << "Test 4: basic geometric series basis test" << endl;

  controller cont1("test4");
  cout << "permutation_number = " << cont1.permutation_number << endl;
  cout << "permutation_sign = " << cont1.permutation_sign << endl;

  three_body_system tbs(&cont1, 1.0, 1.0, 1.0, 1.0, 0.01, 0.01, 1E8, 1E8);
  cout << "V0s = " << tbs.V0_like << " " << tbs.V0_unlike << endl;
  tbs.generate_states_geometric();
  tbs.regenerate_matrices();
  tbs.diagonalise();

  cout << "State info: " << endl;
  cout << tbs.states.back() << endl;

  cout << "Trapping: " << tbs.trapping[0][0] << " " << tbs.trapping[1][1] << " " << tbs.trapping[2][2] << " " << endl;
  cout << "Kinetic: " << tbs.kinetic[0][0] << " " << tbs.kinetic[1][1] << " " << tbs.kinetic[2][2] << " " << endl;

  vector<double> res1;
  cont1.ang->print_angular_states();




}


void rng_test(){
  cout << "Performing rng tests..." << endl << endl;

  cout << "Random integers from -5 to 23: " << endl;
  int min_int = 0, max_int = 0;
  for (int i=1; i<=100; i++){
    int new_int = random_int(-5, 23);
    if (new_int > max_int) max_int = new_int;
    if (new_int < min_int) min_int = new_int;
    cout << setw(3) << new_int << " ";
    if (i % 5 == 0) cout << endl;
  }

  cout << "Max = " << setw(3) << max_int << ", min = " << setw(3) << min_int << endl << endl;

  cout << "Random doubles between -1000.666 and 10000.332:" << endl;
  double min_double = 0.0, max_double = 0.0;
  cout.precision(4);
  for (int i=1; i<=100; i++){
    double new_double = random_double(-1000.888, 10000.332);
    if (new_double > max_double) max_double = new_double;
    if (new_double < min_double) min_double = new_double;
    cout << setw(7) << new_double << " ";
    if (i % 5 == 0) cout << endl;
  }
  cout << "Max = " << setw(7) << max_double << ", min = " << setw(7) << min_double << endl << endl;

}

void compare_vectors(vector<double> &v1, vector<double> &v2){

  if (v1.size() != v2.size()){
    cout << "VECTORS NOT SAME SIZE" << endl;
    return;
  }

  for (unsigned int i=0; i<v1.size(); i++){
    if (v1[i] != v2[i]){
      cout << v1[i] << " " << v2[i] << endl;
    }
  }
}

void testP(angular_subsystem *a){
  cout << "P integral test: " << endl;
  vector<double> P_compare;
  cout << read_array_from_file(P_compare, "P_compare", "hi", a->cache_P.size()) << endl;
  cout << "P_compare.size() = " << P_compare.size() << endl;
  cout << "P.size() = " << a->cache_P.size() << endl;

  compare_vectors(P_compare, a->cache_P);
}

void test7(){
  controller cont1("test5");
  cout << "permutation_number = " << cont1.permutation_number << endl;
  cout << "permutation_sign = " << cont1.permutation_sign << endl;
  three_body_system tbs(&cont1, 1.0, 1.0, 2.0, 1.0, 0.01, 0.01, 1E-8, 1E-8);
  tbs.focus_states_near_r0 = true;
  cont1.ang->print_angular_states();

  while (true){
    tbs.states.push_back(tbs.generate_random_state(0));
    cout << "new state " << tbs.states.back().ang->index << endl;
    tbs.regenerate_matrices();
    tbs.diagonalise();

  }
  //tbs.stochastic_loop();
}

void hyper_test(string name){

  controller c(name);

  three_body_system tbs(&c, 1.0, 1.0 , 1.0, 1.0, 0.01, 0.01, 1E8, 1E8);

  int i1 = 0;
  int i2 = 5;
  tbs.generate_states_geometric();
  tbs.regenerate_matrices();
  tbs.diagonalise();
  cout << "States: " << endl;
  cout << tbs.states[i1];
  cout << tbs.states[i2];
  cout << "E = " << c.ang->E(1,1,0,1,1,0,0,0,0) << endl;

  double rmin = 0.0;
  double rmax = 0.006;
  int nr = 200;
  double rgrad = (rmax - rmin) / double(nr);

  ofstream op("htest");
  for (int i = 0; i < nr; i++){
    double r = rmin + rgrad * i;
    op << r << " " << tbs.calc_sf_channel_element(r, tbs.states[i1], tbs.states[i2], 0) / double(c.permutation_number) << endl;
  }
  op.close();

  vector<double> tar;
  tbs.calc_matrix_elements(tar, tbs.states[i1], tbs.states[i2]);
  cout << "S = " << tar[0] << endl;

  gsl_function F;
  F.function = &(hyper_helper_function);
  vector<double> transformed_state;
  tbs.calc_transformed_state(transformed_state, tbs.states[i1], tbs.states[i2], 3);
  double eta = tbs.length_hyperspherical;

  double h = 0.102;
  double rho_max = h / sqrt(eta);

  hyper_integrand_helper hih(&transformed_state, h, eta, 0, 0, 0, 0, 0, 1.0);
  F.params = &hih;

  //h = 5.7;
  rmin = 0.0;
  rmax = rho_max;
  rmax = 0.003;
  nr = 200;
  double pf = 4.0 * Pi / tbs.states[i1].norm / tbs.states[i2].norm * sqrt(Pi / 2.0) * tbs.ang->E(0,0,0,0,0,0,0,0,0) * tbs.gaussian_integral(tbs.states[i1].alpha_com + tbs.states[i2].alpha_com, 0);
  rgrad = (rmax - rmin) / double(nr);
  cout << "rho_max = " << rho_max << endl;

  ofstream op2("htest2");
  for (int i = 0; i < nr; i++){
    double r = rmin + rgrad * i;
    op2 << r << " " << hyper_helper_function(r, (void *) &hih) << endl;
  }
  cout << "summand = " << tbs.calc_sf_channel_contribution(h, tbs.states[i1], tbs.states[i2], 0, 0, 0, 1) << endl;

  cout << "integral = " << integrate(F, 0.0, rmax, 1E-18, 1E-18) * pf << endl;
  cout << "pf = " << pf << endl;

}
