#include "three_body_system.h"
#include "controller.h"
#include "angular_subsystem.h"
#include "findV0.h"

three_body_system::three_body_system(controller *cont_, double m1_, double m2_, double omega1_, double omega2_, double r0_like_, double r0_unlike_, double a0_like_, double a0_unlike_) : name(cont_->name), cont(cont_), ang(cont_->ang), m1(m1_), m2(m2_), omega1(omega1_), omega2(omega2_), r0_like(r0_like_), r0_unlike(r0_unlike_), a0_like(a0_like_), a0_unlike(a0_unlike_), lim1(cont_->ang->lim1), lim_tr(cont_->ang->lim_tr), lim_kappa(cont_->ang->lim_kappa){




  initialise();
}

void three_body_system::initialise(){
  bb = 0;
  cc = 3;
  only_2body = 0;
  masses[0] = m1;
  hyper_rescale_param = 15.0;
  masses[1] = m2;
  masses[2] = m1;
  omegas[0] = omega1;
  omegas[1] = omega2;
  omegas[2] = omega1;
  channels_allowed = vector<int>{-1, -2, -3, 1, 2, 3};
  initialise_transformations();
  generate_coordinate_matrices();
  output_coordinate_matrices();
  generate_factorials();
  generate_B_coef_caches();

  double prefactor_like = 1.0 / m1;
  double prefactor_unlike = 0.5 * (m1+m2) / m1 / m2;
  r0_like *= sqrt(prefactor_like);
  r0_unlike *= sqrt(prefactor_unlike);
  V0_like = prefactor_like * calc_V0(r0_like, a0_like);
  V0_unlike = prefactor_unlike * calc_V0(r0_unlike, a0_unlike);

  cout << "V0_unlike = " << V0_unlike << endl;
  cout << "masses = " << m1 << " " << m2 << endl;

  basis_max_r = cont->basis_max_r;
  basis_max_rho = cont->basis_max_rho;
  basis_max_com = cont->basis_max_com;
  double r0_minimal = (cont->permutation_sign == 1) ? min(r0_like, r0_unlike) : r0_unlike;
  if (!cont->use_interaction) r0_minimal = 1.0;
  basis_min_r_scaled = (cont->scale_min_to_r0_r) ? cont->basis_min_r * r0_minimal : cont->basis_min_r;
  basis_min_rho_scaled = (cont->scale_min_to_r0_rho) ? cont->basis_min_rho * r0_minimal : cont->basis_min_rho;
  basis_min_com_scaled = (cont->scale_min_to_r0_com) ? cont->basis_min_com* r0_minimal : cont->basis_min_com;

  length_hyperspherical = pow(m1 + m2, 2) / m2 / (m2 + 2.0 * m1);
  report_length_scales();
}

void three_body_system::generate_factorials(){
  int limit = 3*ang->lim10p3;

  gamma.reserve(limit);
  gamma.push_back(tgamma(1.5));
  for (int n=1; n<limit; n++) gamma.push_back(gamma[n-1]*(n+0.5));

  gamma2.reserve(limit);
  for (int n=0; n<limit; n++) gamma2.push_back(tgamma((3.0+n)/2.0));

  fact.reserve(limit);
  fact.push_back(1.0);
  for (int n=1; n<limit; n++) fact.push_back(fact[n-1]*n);

  fact2.reserve(limit);
  fact2.push_back(1.0);
  fact2.push_back(1.0);
  for (int n=2; n<limit; n++) fact2.push_back(fact2[n-2]*n);
}

ostream& operator<<(ostream &os, three_body_system const &t){
  os << t.m1 << " " << t.m2 << " | " << t.omega1 << " " << t.omega2 << " | " << t.r0_like << " " << t.r0_unlike << endl;
  return os;

}

void three_body_system::reset(){
  states = vector<state_full>();


}
