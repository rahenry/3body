#include "radial_subsystem.h"
#include "angular_subsystem.h"
#include "states.h"
#include "three_body_system.h"


radial_subsystem::radial_subsystem(three_body_system *tbs_, angular_subsystem *ang_, double m1_, double m2_, double omega1_, double omega2_, double r0_like_, double r0_unlike_) : tbs(tbs_), ang(ang_), m1(m1_), m2(m2_), omega1(omega1_), omega2(omega2_), r0_like(r0_like_), r0_unlike(r0_unlike_) {
  /*scale_min_r = tbs->basis_min_r * pow(trapping[0][0] / kinetic[0][0], -0.25);
  if (scale_min_to_r0_r) scale_min_r *= r0;
  scale_max_r = tbs->basis_max_r * pow(trapping[0][0] / kinetic[0][0], -0.25);
  scale_min_rho = tbs->basis_min_rho * pow(trapping[1][1] / kinetic[1][1], -0.25);
  if (scale_min_to_r0_rho) scale_min_rho *= r0;
  scale_max_rho = tbs->basis_max_rho * pow(trapping[1][1] / kinetic[1][1], -0.25);
  scale_min_com = tbs->basis_min_com * pow(trapping[2][2] / kinetic[2][2], -0.25);
  if (scale_min_to_r0_com) scale_min_com *= r0;
  scale_max_com = tbs->basis_max_com * pow(trapping[2][2] / kinetic[2][2], -0.25);*/

}


/*state_full radial_subsystem::generate_new_state(){


  new_state.alpha_r = rng_alpha_r(rand_engine);
  new_state.alpha_rho = rng_alpha_rho(rand_engine);
  new_state.alpha_R = rng_alpha_R(rand_engine);


  //return state_full(...);




}*/
