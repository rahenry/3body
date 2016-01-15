#include "three_body_system.h"
#include "angular_subsystem.h"
#include "controller.h"
#include <gsl/gsl_sf_laguerre.h>
#include "integration.h"
#include "structure_factors.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_min.h>

void three_body_system::print_structure_factors(){

  if (cont->sf_channels.size() == 0) return;

  int sf_index = 0;
  for (auto & ch : cont->sf_channels){
    cout << "Calculating sf " << ch << endl;
    double oscillator_length = 1.0;
    //if (ch != 0){
    if (0){
      int cind = channel_index(ch);
      oscillator_length = pow(trapping[cind][0][0] / kinetic[cind][0][0], -0.25) / sqrt(2.0);
    }
    //else oscillator_length = length_hyperspherical;
    else oscillator_length = pow(trapping[channel_index(3)][0][0] / kinetic[channel_index(3)][0][0], -0.25) / sqrt(2.0);
    ofstream sf_output(name + ".sf" + to_string(ch));
    for (unsigned int i = 0; i < sf_radii.size(); i++){
      //double temp1 = (ch == 0) ? oscillator_length : 1.0;
      double temp1 = 1.0;
      sf_output << temp1 * sf_radii[i] << " " << sf_results[sf_index][i] * oscillator_length / double(cont->permutation_number) << endl;
    }
    cout << endl;
    sf_index++;


  }

}


void three_body_system::calc_structure_factors_channels(){

  if (cont->sf_channels.size() == 0) return;
  cout << "Calculating structure factors: ";
  for (auto & ch : cont->sf_channels){
    sf_results.push_back(vector<double>());
    cout << ch << " ";
  }
  cout << endl;

  double sf_grad = (cont->sf_maximum - cont->sf_minimum) / double(cont->sf_n_points);
  for (int sf_point_index = 0; sf_point_index < cont->sf_n_points; sf_point_index++){
    cout << double(sf_point_index) / double(cont->sf_n_points) << " " << flush;
    double rP = sf_point_index * sf_grad + cont->sf_minimum;
    sf_radii.push_back(rP);

    int sf_ch_index = 0;
    for (auto & ch : cont->sf_channels){
      double oscillator_length = 1.0;
      if (ch != 0){
	int cind = channel_index(ch);
	oscillator_length = pow(trapping[cind][0][0] / kinetic[cind][0][0], -0.25) / sqrt(2.0);
      }
      //else oscillator_length = length_hyperspherical;
      else oscillator_length = pow(trapping[channel_index(3)][0][0] / kinetic[channel_index(3)][0][0], -0.25);
      //else oscillator_length = 1.0;
      vector<double> sf_matrix_elements;
      for (auto & s1 : states){
        for (auto & s2 : states){
          sf_matrix_elements.push_back(
              calc_sf_channel_element(oscillator_length*rP, s1, s2, ch)
              );

        }
      }

      double res = 0.0;
      for (int i=0; i<basis_size; i++){
        for (int j=0; j<basis_size; j++){
          res += evecs[i] * evecs[j] * sf_matrix_elements[i + basis_size * j];
        }
      }

      sf_results[sf_ch_index++].push_back(res);

      if (test_structure_factors){
	sf_test_matrices.push_back(vector<double>());
	for (int i=0; i<basis_size; i++){
	  for (int j=0; j<basis_size; j++){
	    sf_test_matrices.back().push_back(sf_matrix_elements[i + basis_size * j]);
	  }
	}
      }
    }
  }
  cout << "... " << endl;
  compute_sf_test_results();
  print_structure_factors();

}

void three_body_system::compute_sf_test_results(){
  if (!test_structure_factors) return;



}

double three_body_system::calc_sf_channel_element(double rP, state_full s1, state_full s2, int channel_working){
  double res = 0.0;

  res += calc_sf_channel_loop_transformations(rP, s1, s2, channel_working, 1.0);

  if (cont->permutation_number == 2){
    s1.permute();
    res += calc_sf_channel_loop_transformations(rP, s1, s2, channel_working, cont->permutation_sign);

    //if (!(abs(channel_working) == 5)){
    if (channel_working != 0){
      s2.permute();
      res += calc_sf_channel_loop_transformations(rP, s1, s2, channel_working, 1.0);

      s1.permute();
      res += calc_sf_channel_loop_transformations(rP, s1, s2, channel_working, cont->permutation_sign);
    }
  }
  
  return res;



}

double three_body_system::calc_sf_channel_loop_transformations(double rP, state_full &s1, state_full &s2, int channel_working, double prefactor){
  double res = 0.0;
  bool quiet = 1;

  int c = (channel_working == 0) ? 3 : channel_working;

  for (unsigned int transformation_index1 = 0; transformation_index1 < s1.ang->allowed_transformations.size(); transformation_index1++){
    if (!quiet) cout << "B1 = B1(" << transformation_index1 << ", " << channel_working << ")" << endl << flush;
    double B1 = B(s1, transformation_index1, c);

    if (!quiet) cout << B1 << endl << flush;
    if (B1 == 0.0) continue;

    for (unsigned int transformation_index2 = 0; transformation_index2 < s2.ang->allowed_transformations.size(); transformation_index2++){
      double B2 = B(s2, transformation_index2, c);
      if (B2 == 0.0) continue;

      res += calc_sf_channel_contribution(rP, s1, s2, transformation_index1, transformation_index2, channel_working, prefactor * B1 * B2);

    }
  }


  return res;


}


double three_body_system::calc_sf_channel_contribution(double rP, state_full &s1, state_full &s2, int transformation_index1, int transformation_index2, int channel_working, double prefactor){
  double res = 0.0;
  int c = (channel_working == 0) ? 3 : channel_working;

  state_angular *a1 = &ang->states_transformed[s1.ang->allowed_transformations[transformation_index1]];
  state_angular *a2 = &ang->states_transformed[s2.ang->allowed_transformations[transformation_index2]];

  if (a1->lcom != a2->lcom) return 0.0;

  int cind = channel_index(c);

  prefactor *= 4.0 * Pi / s1.norm / s2.norm * sqrt(Pi / 2.0);

  for (int kappa=0; kappa < get_kappa_limit(a1, a2); kappa++){
    double E = ang->E(kappa, kappa, 0, a2->lr, a2->lrho, a2->lcoup, a1->lr, a1->lrho, a1->lcoup);
    if (E == 0.0) continue;

    vector<double> transformed_state;
    calc_transformed_state(transformed_state, s1, s2, c);

    double nr_double = a1->nr + a2->nr + (a1->lr + a2->lr - kappa) / 2.0 + 0.01;
    double nrho_double  = a1->nrho + a2->nrho + (a1->lrho + a2->lrho - kappa) / 2.0 + 0.01;
    int nr = floor(nr_double);
    int nrho = floor(nrho_double);
    double com0 = gaussian_integral(s1.alpha_com + s2.alpha_com, a1->lcom + a2->lcom);

    if (channel_working != 0) res += com0 * E * calc_sf_channel_summand(kappa, nr, nrho, transformed_state, rP);
    else res += com0 * E * calc_sf_hyper_summand(kappa, nr, a1->lr, nrho, a1->lrho, transformed_state, rP);

  }

  res *= prefactor;

  return res;



}

double three_body_system::calc_sf_hyper_summand(int kappa, int nr, int lr, int nrho, int lrho, vector<double> transformed_state, double h){
  gsl_function F;
  F.function = &(hyper_helper_function);
  double eta = length_hyperspherical;
  double rho_max = h / sqrt(eta);
  hyper_integrand_helper hih(&transformed_state, h, eta, kappa, lr, nr, lrho, nrho, hyper_rescale_param);
  F.params = &hih;

  double m = 0.99*pow(rho_max, hyper_rescale_param);
  F.function = &(hyper_helper_function_scaled);
  m = minimise(F, 0.0, pow(rho_max, hyper_rescale_param), m, 1E-12, 1E-4);

  double mm = pow(m, 1.0 / hyper_rescale_param);
  double lim = min(100 * mm, rho_max);
  //cout << "rho_max = " << rho_max << ", min = " << mm << ", lim = " << lim << endl;


  F.function = &(hyper_helper_function);
  double res = integrate(F, 0.0, lim, 0.0, 1E-6);
  //double res = -integrate_to_inf(F, 0.0, 1E-8, 1E-6);


  return 2.0 * res;

}

double hyper_helper_function_scaled(double rho, void * params){
  hyper_integrand_helper *helper = (hyper_integrand_helper *) params;
  return -abs(hyper_helper_function(pow(rho, 1.0 / helper->rescale_param), params));
}

double hyper_helper_function(double rho, void * params){
  hyper_integrand_helper *helper = (hyper_integrand_helper *) params;
  double u = (*helper->transformed_state)[0];
  double v = (*helper->transformed_state)[2];
  double w = (*helper->transformed_state)[1];
  double h = helper->h;

  double rsqr = h * h - (helper->eta) * rho * rho;
  //cout << h << " " << rsqr << " " << rho << endl;
  if (rsqr < 0){
    //cout << h << " " << helper->eta << " " << rho << endl;
    return 0.0; // in this case there is no point that satisfies the delta function
  }


  //double r = sqrt(abs(rsqr));
  double r = sqrt(rsqr);
  //cout << rsqr << endl;
  double arg = rho * r * abs(w);
 
  double temp1 = -0.5 * (u * r * r + v * rho * rho);
  double temp2 = pow(r, 2*helper->nr + helper->kappa + 1) * pow(rho, 2*helper->nrho + helper->kappa + 2) * h;
  double res = 0.0;
  //cout << temp1 << " " << temp2 << " " << arg << endl;
  if (abs(w/u) < 1E-8){ 
    if (helper->kappa != 0){ 
      return 0.0;
    }
    res = exp(temp1) * temp2 / sqrt(Pi / 2.0);
  }
  else{
    gsl_sf_result bes1;

    int status = gsl_sf_bessel_il_scaled_e(helper->kappa, arg, &bes1);
    if (status != 0){
      bes1.val = 1.0 / arg;
    }
    res = bes1.val * exp(abs(arg) + temp1) * temp2 * epsilon(helper->kappa, w) * sqrt(2.0 * helper->kappa + 1.0) / sqrt(Pi / 2.0);
  }
  //cout << res << " " ;
  return res;
}


double three_body_system::calc_sf_channel_summand(int kappa, int nr, int nrho, vector<double> transformed_state, double rP){
  double u = transformed_state[0];
  double v = transformed_state[2];
  double w = transformed_state[1];

  return
    epsilon(kappa, w) * sqrt(2.0*kappa + 1.0) * fact2[2*nrho] * pow(abs(w), kappa)
    / pow(v, 1.5 + nrho + kappa) * pow(rP, 2 + 2*nr + 2*kappa) * exp(-0.5 * (u - w*w/v) * rP * rP) 
    * gsl_sf_laguerre_n(nrho, double(kappa) + 0.5, -rP * rP * w * w / v / 2.0);


}


