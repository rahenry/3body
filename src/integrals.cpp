#include "three_body_system.h"
#include "angular_subsystem.h"
#include "states.h"
#include "gsl/gsl_blas.h"
#include "controller.h"
#include "utility.h"

void three_body_system::calc_matrix_element_contributions(vector<double> &tar, state_full &s1, state_full &s2, int transformation_index1, int transformation_index2, double prefactor, int channel_working, bool calc_interaction_only){
  bool quiet = 1;
  state_angular *a1 = &ang->states_transformed[s1.ang->allowed_transformations[transformation_index1]];
  if (channel_working == s1.channel){
    a1 = s1.ang;
  }
  state_angular *a2 = &ang->states_transformed[s2.ang->allowed_transformations[transformation_index2]];
  double overlap = 0.0;
  double H_trapping = 0.0;
  double H_kinetic = 0.0;
  double H_interaction = 0.0;

  int cind = channel_index(channel_working);

  prefactor *= 4.0 * Pi / s1.norm / s2.norm;
  //if (a1->index == 4 && a2->index==2 && !calc_interaction_only) cout << get_kappa_limit(a1, a2) << " " << calc_interaction_only << endl;
        //cout << a1->index << " " << a2->index << endl;
  if (!quiet){
    cout << "Angular states: " << endl << *a1 << endl << *a2 << endl;
    cout << "kappa_max = " << get_kappa_limit(a1, a2) << endl;
  }
  for (int kappa=0; kappa < get_kappa_limit(a1, a2); kappa++){
    if (kappa > ang->lim_kappa) cout << "ALKJDSLJKSADLKJADS";
    double E = ang->E(kappa, kappa, 0, a2->lr, a2->lrho, a2->lcoup, a1->lr, a1->lrho, a1->lcoup);
    if (!quiet){
      cout << kappa << ": " << E << endl;
    }
    //if (!cont->use_trapping_imbal && E == 0.0) continue;

    vector<double> transformed_state;
    calc_transformed_state(transformed_state, s1, s2, channel_working);

    if (!quiet){
      cout << "Transformed sums: " << endl;
      print_vector(transformed_state);
    }

    double nr_double = a1->nr + a2->nr + (a1->lr + a2->lr - kappa) / 2.0 + 0.01;
    double nrho_double  = a1->nrho + a2->nrho + (a1->lrho + a2->lrho - kappa) / 2.0 + 0.01;
    int nr = floor(nr_double);
    int nrho = floor(nrho_double);

    if (nr_double < -0.9 || nrho_double < -0.9){
      //cout << "WHAT?" << endl;
      continue;
    }

    double com0 = gaussian_integral(s1.alpha_com + s2.alpha_com, a1->lcom + a2->lcom);
    //cout << calc_interaction_only << " " << cont->use_trapping << " " << cont->use_trapping_imbal << endl;
    if (!calc_interaction_only && cont->use_trapping && cont->use_trapping_imbal){
      // we use s1.ang since these will only be calculated in cases where we transform s2 to s1's channel
  //if (a1->index == 4 && a2->index==2) cout << "ASDAS" << get_kappa_limit(a1, a2) << endl;
      double F = ang->F(a1, a2, kappa);
      double H_r = ang->H_r(a1, a2, kappa);
      double H_rho = ang->H_rho(a1, a2, kappa);
      double com_imb = gaussian_integral(s1.alpha_com + s2.alpha_com, a1->lcom + a2->lcom + 1);
      //cout << "F = " << F << endl;
      //if (H_r != 0.0) cout << "H_r = " << H_r << endl;
      //cout << "H_rho = " << H_rho << endl;

      if (F != 0.0){
        //cout << "F: " << nr_double << " " << nrho_double << endl;
        H_trapping += trapping[cind][1][0] * prefactor * F * com0 * 
          bessel_integral(transformed_state, kappa, int(nr_double + 0.5), int(nrho_double + 0.5), 0.0);
      }
      if (H_r != 0.0){
        //cout << "H_r: " << nr_double << " " << nrho_double << endl;
        H_trapping += trapping[cind][0][2] * prefactor * H_r * com_imb * 
          bessel_integral(transformed_state, kappa, int(nr_double + 0.5), nrho, 0.0);
      }
      if (H_rho != 0.0){ 
        //cout << "H_rho: " << nr_double << " " << nrho_double << endl;
        H_trapping += trapping[cind][1][2] * prefactor * H_rho * com_imb * 
          bessel_integral(transformed_state, kappa, nr, int(nrho_double + 0.5), 0.0);
      }
    }

    vector<double> angular_couplings;
    calc_angular_couplings(angular_couplings, a1, a2);
    double cg2 = angular_couplings[1];
    E *= cg2;
    if (E == 0.0) continue;

    double bessel00 = bessel_integral(transformed_state, kappa, nr, nrho, 0.0);
    double bessel10 = bessel_integral(transformed_state, kappa, nr+1, nrho, 0.0);
    double bessel01 = bessel_integral(transformed_state, kappa, nr, nrho+1, 0.0);
    double com1 = gaussian_integral(s1.alpha_com + s2.alpha_com, a1->lcom + a2->lcom + 2);

    if (calc_interaction_only){
      double V0 = (abs(channel_working)==2) ? V0_like : V0_unlike;
      double r0 = (abs(channel_working)==2) ? r0_like : r0_unlike;
      H_interaction += V0 * prefactor * E * com0 * 
	bessel_integral(transformed_state, kappa, nr, nrho, pow(r0, -2));
      continue;
    }

    overlap += prefactor * E * bessel00 * com0;

    if (cont->use_trapping){
      H_trapping += 0.5 * prefactor * E * (
          trapping[cind][0][0] * bessel10 * com0 +
          trapping[cind][1][1] * bessel01 * com0 +
          trapping[cind][2][2] * bessel00 * com1);
    }
    if (cont->use_kinetic){
      double temp1 = kinetic[cind][0][0] * s1.alpha_r * (3.0+2.0*a1->lr) +
        kinetic[cind][1][1] * s1.alpha_rho * (3.0+2.0*a1->lrho) + 
        kinetic[cind][2][2] * s1.alpha_com * (3.0+2.0*a1->lcom);
      H_kinetic += -0.5 * prefactor * E * (
          kinetic[cind][0][0] * pow(s1.alpha_r, 2) * bessel10 * com0 +
          kinetic[cind][1][1] * pow(s1.alpha_rho, 2) * bessel01 * com0 +
          kinetic[cind][2][2] * pow(s1.alpha_com, 2) * bessel00 * com1 -
          temp1 * bessel00 * com0);
    }

  }

  if (!quiet) cout << "Total overlap contribution = " << overlap << endl << "_______" << endl;
  tar[0] += overlap;
  tar[1] += H_trapping;
  tar[2] += H_kinetic;
  tar[3] += H_interaction;

}

double three_body_system::bessel_integral(vector<double> &transformed_state, int kappa, int nr, int nrho, double interaction){
  double w = transformed_state[1];
  if (w == 0.0 && kappa != 0) return 0.0;

  double u = transformed_state[0] + interaction;
  double v = transformed_state[2];

  double w22v = w*w/2.0/v;
  double umw2v = 0.5 * (u - w*w/v);

  double sum1 = 0.0;

  for (int k=0; k<=nrho; k++){
    double pow1 = (w22v == 0.0 && k == 0) ? 1.0 : pow(w22v, k);
    sum1 += gamma[k + nr + kappa] / fact[k] / fact[nrho-k] / gamma[k+kappa] * pow1 / pow(umw2v, 1.5+k+nr+kappa);
  }

  double pow2 = (w == 0.0 && kappa == 0.0) ? 1.0 : pow(abs(w), kappa);
  return sqrt(2.0*kappa+1.0) * epsilon(kappa, w) * sqrt(Pi / 8.0) * fact2[2*nrho] * gamma[nrho+kappa] * pow2 * pow(v, -1.5 - nrho - kappa) * sum1;
}

void three_body_system::loop_transformations_H(vector<double> &tar, state_full &s1, state_full &s2, int channel_working, double prefactor, bool calc_interaction_only){

  int quiet = 1;
  for (unsigned int transformation_index1 = 0; transformation_index1 < s1.ang->allowed_transformations.size(); transformation_index1++){
    double B1 = B(s1, transformation_index1, channel_working);
    if (!quiet) cout << transformation_index1 << ": B1 = " << B1 << endl;
    if (B1 == 0.0) continue;

    for (unsigned int transformation_index2 = 0; transformation_index2 < s2.ang->allowed_transformations.size(); transformation_index2++){
      double B2 = B(s2, transformation_index2, channel_working);
      if (!quiet) cout << transformation_index2 << ": B2 = " << B2 << endl;
      if (B2 == 0.0) continue;
      
      calc_matrix_element_contributions(tar, s1, s2, transformation_index1, transformation_index2, prefactor * B1 * B2, channel_working, calc_interaction_only);
    }
  }
}

void three_body_system::calc_matrix_elements(vector<double> &tar, state_full s1, state_full s2){
  // States are passed by value so they can be permuted in-place
  bool quiet = 1;
  // tar[0] = overlap
  // tar[1] = kinetic
  // tar[2] = trapping
  // tar[3] = interaction
  tar = vector<double>(4, 0.0);

  //loop_transformations_H(tar, s1, s2, s1.channel, 1.0, false);
  loop_transformations_H(tar, s1, s2, s1.channel, 1.0, false);

  if (cont->use_interaction){
    loop_transformations_H(tar, s1, s2, 3, 1.0, true);
    if (!only_2body) loop_transformations_H(tar, s1, s2, -1, 1.0, true);
    //if (cont->permutation_sign == 1) loop_transformations_H(tar, s1, s2, 2, 1.0, true);
  }

  if (cont->permutation_number > 1 && !only_2body){
    s1.permute();
    loop_transformations_H(tar, s1, s2, s1.channel, cont->permutation_sign, false);
    if (cont->use_interaction){
      loop_transformations_H(tar, s1, s2, 3, cont->permutation_sign, true);
      loop_transformations_H(tar, s1, s2, -1, cont->permutation_sign, true);
      if (cont->permutation_sign == 1) loop_transformations_H(tar, s1, s2, 2, cont->permutation_sign, true);
    }
  }
}

int three_body_system::get_kappa_limit(state_angular *a1, state_angular *a2){
  int nr = 2*(a1->nr + a2->nr) + a1->lr + a2->lr;
  int nrho = 2*(a1->nrho + a2->nrho) + a1->lrho + a2->lrho;
  return min(nr, nrho) + 2;
}

double three_body_system::gaussian_integral(double a, int b){
  return 0.5 * gamma2[b] / pow(a * 0.5, 0.5*b + 1.5);
}

double gaussian_integral(double a, int b){
  return 0.5 * tgamma(1.5+0.5*b) / pow(a * 0.5, 0.5*b + 1.5);
}

double epsilon(int kappa, double w){
  if (w > 0.0) return 1.0;
  else return pow(-1.0, kappa);
}

void three_body_system::calc_angular_couplings(vector<double> &tar, state_angular *a1, state_angular *a2){
  tar = vector<double>();
  double cg1 = 0.0, cg2 = 0.0;

  bool match_all = false, match_com = false, match_r = false, match_rho = false;
  if (a2->lr == a1->lr) match_r = true;
  if (a2->lrho == a1->lrho) match_rho = true;
  if (a2->lcom == a1->lcom){
    match_com = true;
    if (match_r && match_rho){
      match_all = true;
    }
  }

  //double trapImbInt1 = 0.0, trapImbInt2 = 0.0, trapImbInt3 = 0.0;
  double sum1 = 0.0, sum2 = 0.0;
  for (int mcom2=-a2->lcom; mcom2<=a2->lcom; mcom2++){
    for (int mcom1=-a1->lcom; mcom1<=a1->lcom; mcom1++){
      int mcoup1 = -mcom1;
      int mcoup2 = -mcom2;
      if (abs(mcoup1) > a1->lcoup || abs(mcoup2) > a2->lcoup) continue;
      double cg1 = ang->cg(a1->lcoup, mcoup1, a1->lcom, mcom1, ang->ltotal, 0);
      double cg2 = ang->cg(a2->lcoup, mcoup2, a2->lcom, mcom2, ang->ltotal, 0);

      if (match_com && mcom1==mcom2){
	sum2 += cg1*cg2;
      }
      for (int mrho1=-a1->lrho; mrho1<=a1->lrho; mrho1++){
	for (int mrho2=-a2->lrho; mrho2<=a2->lrho; mrho2++){
	  int mr1 = mcoup1 - mrho1;
	  int mr2 = mcoup2 - mrho2;
	  if (abs(mr1) > a1->lr || abs(mr2) > a2->lr) continue;
	  double cg3 = ang->cg(a1->lr, mr1, a1->lrho, mrho1, a1->lcoup, mcoup1);
	  double cg4 = ang->cg(a2->lr, mr2, a2->lrho, mrho2, a2->lcoup, mcoup2);

	  double cgprod = cg1*cg2*cg3*cg4;

	  if (match_all && (mcom1==mcom2 && mrho1==mrho2)){
	    sum1 += cgprod;
	  }

	  /*if (useTrappingImbal){
	    if (match_com && mcom==mcomP) trapImbInt1 += cgprod * AMIntegral_P(lrP,mrP,lrhoP,mrhoP,lr,mr,lrho,mrho);
	    if (match_r && mr==mrP) trapImbInt2 += cgprod * AMIntegral_P(lrhoP,mrhoP,lcomP,mcomP,lrho,mrho,lcom,mcom);
	    if (match_rho && mrho==mrhoP) trapImbInt3 += cgprod * AMIntegral_P(lrP,mrP,lcomP,mcomP,lr,mr,lcom,mcom);
	  }*/
	}
      }
    }
  }

  tar.push_back(sum1);
  tar.push_back(sum2);



}


