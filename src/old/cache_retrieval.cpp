#include "angular_subsystem.h"
#include "states.h"

double angular_subsystem::cg(int l1, int m1, int l2, int m2, int l3, int m3){
  return cache_cg[
    (m3 + m_offset_cg) + lim2_cg * (
      l3 + lim_cg * (
	(m2 + m_offset_cg) + lim2_cg * (
	  l2 + lim_cg * (
	    (m1 + m_offset_cg) + lim2_cg * (
	      l1
	      )))))];
}

double angular_subsystem::P(int l1, int m1, int l2, int m2, int l3, int m3, int l4, int m4){

  if (l1+l3==0 || l2+l4==0) return 0.0;
  if (abs(l1-l3) > 1 || abs(l2-l4) > 1) return 0.0;
  if (m1 - m3 + m2 - m4 != 0) return 0.0;
  if (abs(m1-m3) > 1) return 0.0;

  return cache_P[
    (m4-m2+1) + 3 * (
      (l4-l2+1) + 3 * (
        (m3-m1+1) + 3 * (
          (l3-l1+1) + 3 * (
            (m2 +lim_special-1) + (2*lim_special-1) * (
              l2 + lim_special * (
                (m1 +lim_special-1) + (2*lim_special-1) * (
                  l1
		  )))))))];
}


double angular_subsystem::Q(state_angular &s1, state_angular &s2, int a, int b, int mcom){
  return cache_Q[
    (mcom + lmax) + lim2 * (
      b + lim_special * (
	a + lim_special * (
	  s2.index + n_states_transformed *
	    s1.index
	  )))];
}

double angular_subsystem::F(state_angular &s1, state_angular &s2, int kappa){
  return cache_F[
    kappa + lim_kappa * (
	s2.index + n_states_transformed * s1.index)
    ];
}

double angular_subsystem::G_r(state_angular &s1, state_angular &s2, int a){
  return cache_G_r[
    a + lim_special * (
	s2.index + n_states_transformed * s1.index)
    ];
}


double angular_subsystem::G_rho(state_angular &s1, state_angular &s2, int b){
  return cache_G_rho[
    b + lim_special * (
	s2.index + n_states_transformed * s1.index)
    ];
}

double angular_subsystem::H_r(state_angular &s1, state_angular &s2, int kappa){
  return cache_H_r[
    kappa + lim_special * (
	s2.index + n_states_transformed * s1.index)
    ];
}

double angular_subsystem::H_rho(state_angular &s1, state_angular &s2, int kappa){
  return cache_H_rho[
    kappa + lim_special * (
	s2.index + n_states_transformed * s1.index)
    ];
}

double angular_subsystem::E(int l1, int l2, int l3, int l4, int l5, int l6, int l7, int l8, int l9){
  return cache_E[
    l9 + lim1 * (
      l8 + lim_special * (
        l7 + lim_special * (
          l6 + lim1 * (
            l5 + lim2 * (
              l4 + lim2 * (
                l3 + lim1 * (
                  l2 + lim_kappa *
                    l1
                  )))))))];
}
          
