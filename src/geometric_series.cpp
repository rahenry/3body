#include "three_body_system.h"
#include "controller.h"
#include "angular_subsystem.h"
#include "utility.h"

void three_body_system::geometric_series(vector<double> &tar, double geo_min, double geo_max, double N){
  tar = vector<double>();
  tar.resize(N);
  double current_term = geo_min;
  double ratio = pow(geo_max / geo_min, 1.0 / (N-1.0));
  for (auto & element : tar){
    element = pow(current_term, -2);
    current_term *= ratio;
  }
}


void three_body_system::generate_states_geometric(){
  vector<double> alphas_r, alphas_rho, alphas_com;
  int channel = 3;
  int cind = channel_index(channel);

  double minimum = cont->scale_min_to_r0_r ? cont->basis_min_r * r0_unlike : cont->basis_min_r;
  double maximum = cont->basis_max_r * pow(trapping[cind][0][0] / kinetic[cind][0][0], -0.125);
  cout << cont->basis_max_r << " " << pow(trapping[cind][0][0] / kinetic[cind][0][0], -0.125) << endl;
  geometric_series(alphas_r, minimum, maximum, 2);
  minimum = cont->scale_min_to_r0_rho ? cont->basis_min_rho * r0_unlike : cont->basis_min_rho;
  maximum = cont->basis_max_rho * pow(trapping[cind][1][1] / kinetic[cind][1][1], -0.125);
  geometric_series(alphas_rho, minimum, maximum, 2);

  if (1 == 1){
    //alphas_com.push_back(pow(trapping[2][2] / kinetic[2][2], 0.25));
    alphas_com.push_back(1.0/kinetic[cind][2][2]);
  }
  else{
    minimum = cont->scale_min_to_r0_com ? cont->basis_min_com * r0_unlike : cont->basis_min_com;
    maximum = cont->basis_max_com * pow(trapping[cind][2][2] / kinetic[cind][2][2], -0.125);
    geometric_series(alphas_com, minimum, maximum, 1);
  }

  states = vector<state_full>();
  int index = 0;
	for (auto & alpha_com : alphas_com){
      for (auto & alpha_rho : alphas_rho){
    for (auto & alpha_r : alphas_r){
  for (auto & angular: ang->states){
          states.push_back(
              state_full(index++, &angular, alpha_r, alpha_rho, alpha_com, 3));
        }
      }
    }
  }
  cout << "Generated geometric series states, total: " << alphas_r.size() << " x " << alphas_rho.size() << " x " << alphas_com.size() << " x " << ang->states.size() << " = " << states.size() << endl;

  if (true){
    print_vector_to_file(alphas_r, name + ".info.alphas_r");
    print_vector_to_file(alphas_rho, name + ".info.alphas_rho");
    print_vector_to_file(alphas_com, name + ".info.alphas_com");
  }

}
