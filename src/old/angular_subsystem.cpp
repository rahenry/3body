#include "angular_subsystem.h"
#include "three_body_system.h"
#include <gsl/gsl_sf_coupling.h>
#include "utility.h"
#include "states.h"

angular_subsystem::angular_subsystem(three_body_system *tbs) :
  name(tbs->name), ltotal(tbs->ltotal), lmax_r(tbs->lmax_r), lmax_rho(tbs->lmax_rho), lmax_coup(tbs->lmax_coup), lmax_com(tbs->lmax_com), single_angular_state(tbs->single_angular_state), use_trapping_imbal(tbs->use_trapping_imbal), parity(tbs->parity)
{
  lmax = max(
      max(lmax_r, lmax_rho), max(lmax_coup, lmax_com)
      );
  lim1 = lmax + 1;
  lim2 = 2 * lmax + 1;
  lim3 = 3 * lmax + 1;
  lim4 = 4 * lmax + 1;
  lim5 = 5 * lmax + 1;
  lim6 = 6 * lmax + 1;
  lim12 = 12 * lmax + 1;
  lim10p3 = 10 * lmax + 4;
  lim_kappa = 3 * lmax + 2;
  lim_special = lim2 + 2;
  lim_cg = max(lim_kappa,lim_special) + 0;
  lim2_cg = 2*lim_cg - 1;
  m_offset_cg = lim_cg - 1;
  lim_coup = lmax_coup + 1;


  generate_states();
  generate_states_transformed();
  generate_caches();
}

void angular_subsystem::generate_states(){
  cout << "Generating angular states... ";

  generate_angular_states(states, lmax_r, lmax_rho, lmax_coup, lmax_com, ltotal, parity, single_angular_state, false, lmax_r+lmax_rho);

  n_states = states.size();

  cout << "total states: " << n_states << endl;
}

void angular_subsystem::generate_states_transformed(){

  cout << "Generating transformed states.. ";

  generate_angular_states(states_transformed, lmax_r*2, lmax_rho*2, lmax_coup, lmax_com, ltotal, parity, single_angular_state, true, lmax_r+lmax_rho);

  n_states_transformed = states_transformed.size();

  for (auto &state : states) state.generate_allowed_transformations(states_transformed);

  cout << "total states: " << n_states_transformed << endl;

}

void angular_subsystem::generate_caches(){
  generate_cg();
  generate_P();
  generate_E();
  generate_Q();
  generate_F();
  generate_G();



}

void angular_subsystem::generate_cg(){

  cout << "Generating CG coefficients... ";
  for (int l1=0; l1<lim_cg; l1++){
    cout << double(l1)/double(lim_cg) << " " << flush;
    for (int m1=-lim_cg+1; m1<lim_cg; m1++){
      for (int l2=0; l2<lim_cg; l2++){
        for (int m2=-lim_cg+1; m2<lim_cg; m2++){
          for (int l3=0; l3<lim_cg; l3++){
            for (int m3=-lim_cg+1; m3<lim_cg; m3++){
              cache_cg.push_back(
                  pow(-1,l1-l2+m3) *
                  sqrt(2.0 * l3 + 1.0) *
                  gsl_sf_coupling_3j(2*l1, 2*l2, 2*l3, 2*m1, 2*m2, -2*m3));
  }}}}}}
  cout << 1.0 << endl;

}

void angular_subsystem::generate_P(){

  int array_length = pow(lim_special, 2) * pow(2 * lim_special + 1, 2) * 81;

  cache_P.reserve(array_length);

  for (int a=0; a<lim_special; a++){
  for (int ma=-lim_special+1; ma<lim_special; ma++){
  for (int b=0; b<lim_special; b++){
  for (int mb=-lim_special+1; mb<lim_special; mb++){
  for (int c=a-1; c<=a+1; c++){
  for (int mc=ma-1; mc<=ma+1; mc++){
  for (int d=b-1; d<=b+1; d++){
  for (int md=mb-1; md<=mb+1; md++){
    if (d >= lim_special || c >= lim_special){
      cache_P.push_back(0.0);
      continue;
    }
    if (c==-1 || d==-1){
      cache_P.push_back(0.0);
      continue;
    }
    if (abs(mc) > c || abs(md) > d){
      cache_P.push_back(0.0);
      continue;
    }

    double cg1 = cg(c,0,1,0,a,0);
    if (cg1 == 0.0){
      cache_P.push_back(0.0);
      continue;
    }
    double cg2 = cg(d,0,1,0,b,0);
    if (cg2 == 0.0){
      cache_P.push_back(0.0);
      continue;
    }

    double res = 0.0;
    for (int m=-1; m<=1; m++){
      double cg3 = cg(a,ma,1,m,c,mc);
      double cg4 = cg(b,mb,1,-m,d,md);
      res += cg3 * cg4 * pow(-1,m);
    }
    cache_P.push_back(res * cg1 * cg2);

  }}}}}}}}

  cout << 1.00 << endl;
}

void angular_subsystem::generate_E(){
  cache_E.reserve(pow(lim1,3) * pow(lim_special,2) * pow(lim2,2) * pow(lim_kappa,2));

  for (int l1=0; l1<lim_kappa; l1++){
  for (int l2=0; l2<lim_kappa; l2++){
  for (int l3=0; l3<lim1; l3++){
  for (int l4=0; l4<lim2; l4++){
  for (int l5=0; l5<lim2; l5++){
  for (int l6=0; l6<lim1; l6++){
  for (int l7=0; l7<lim_special; l7++){
  for (int l8=0; l8<lim_special; l8++){
  for (int l9=0; l9<lim1; l9++){
    double res = gsl_sf_coupling_9j(2*l1, 2*l2, 2*l3, 2*l4, 2*l5, 2*l6, 2*l7, 2*l8, 2*l9);
    res *= sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l4+1.0)*(2.0*l5+1.0)/(2.0*l7+1.0)/(2.0*l8+1.0)) / 4.0 / Pi * cg(l1,0,l4,0,l7,0) * cg(l2,0,l5,0,l8,0);
    cache_E.push_back(res);

  }}}}}}}}}
}

void angular_subsystem::generate_Q(){

  cache_Q.reserve(n_states * n_states_transformed * pow(lim_special, 3));

  for (auto &s1 : states){
    for (auto &s2 : states_transformed){
      for (int a=0; a<lim_special; a++){
	for (int b=0; b<lim_special; b++){
	  for (int mcom=-lim1+1; mcom<lim1; mcom++){
	    double res = 0.0;

	    if (s1.lcom == s2.lcom){
	    for (int mlr2=-s2.lr; mlr2<=s2.lr; mlr2++){
	    for (int mlrho2=-s2.lrho; mlrho2<=s2.lrho; mlrho2++){
	    for (int ma=-a; ma<=a; ma++){
	    for (int mb=-b; mb<=b; mb++){
	      double cg1 = cg(s2.lr, mlr2, s2.lrho, mlrho2, s2.lcoup, -mcom);
	      double cg2 = cg(a, ma, b, mb, s1.lcoup, -mcom);

	      double P1 = P(a, ma, b, mb, s2.lr, mlr2, s2.lrho, mlrho2);

	      res += P1 * cg1 * cg2;
	    }}}}}
	    cache_Q.push_back(res);
	  }
	}
      }
    }
  }
}

void angular_subsystem::generate_F(){

  cache_F.reserve(n_states * n_states_transformed * lim_kappa * pow(lim_special, 2));

  for (auto &s1 : states){
    for (auto &s2 : states_transformed){
      for (int kappa=0; kappa<lim_kappa; kappa++){
	double res = 0.0;

	for (int a=0; a<lim_special; a++){
	  for (int b=0; b<lim_special; b++){
	    double E1 = E(kappa,kappa,0,s2.lr,s2.lrho,s2.lcoup,a,b,s2.lcoup);
	    for (int mcom=-s1.lcom; mcom<=s1.lcom; mcom++){
	      double Q1 = Q(s1,s2,a,b,mcom);
	      res += E1*Q1*cg(s2.lcoup,-mcom,s2.lcom,mcom,ltotal,0)*cg(s1.lcoup,-mcom,s1.lcom,mcom,ltotal,0);
	    }
	  }
	}
	cache_F.push_back(res);
      }
    }
  }
}

void angular_subsystem::generate_G(){
  cout << "Generating \'G\' angular momentum integrals... ";

  cache_G_r.reserve(1);
  cache_G_rho.reserve(1);

  for (auto &s1 : states){
    for (auto &s2 : states_transformed){
      for (int a=0; a<lim_special; a++){
        double res = 0.0;
        for (int mlr2=-s2.lr; mlr2<=s2.lr; mlr2++){
          for (int ma=-a; ma<=a; ma++){
            for (int mcom1=-s1.lcom; mcom1<=s1.lcom; mcom1++){
              for (int mcom2=-s2.lcom; mcom2<=s2.lcom; mcom2++){
                double P1 = P(a, ma, s1.lcom, mcom1, s2.lr, mlr2, s2.lcom, mcom2);
                if (mcom1 + ma != mcom2 + mlr2 ||
                    abs(mcom1 + ma) > s2.lrho ||
                    abs(mcom2 + mlr2) > s2.lrho) continue;

                res += P1 * cg(s2.lr, mlr2, s2.lrho, -mcom2-mlr2, s2.lcoup, -mcom2) *
                  cg(a, ma, s2.lrho, -mcom1-ma, s1.lcoup, -mcom1) *
                  cg(s2.lcoup, -mcom2, s2.lcom, mcom2, ltotal, 0) *
                  cg(s1.lcoup, -mcom1, s1.lcom, mcom1, ltotal, 0);
              }
            }
          }
        }
        cache_G_r.push_back(res);
      }

      for (int b=0; b<lim_special; b++){
        double res = 0.0;
        for (int mlrho2=-s2.lrho; mlrho2<=s2.lrho; mlrho2++){
          for (int mb=-b; mb<=b; mb++){
            for (int mcom1=-s1.lcom; mcom1<=s1.lcom; mcom1++){
              for (int mcom2=-s2.lcom; mcom2<=s2.lcom; mcom2++){
                double P1 = P(b, mb, s1.lcom, mcom1, s2.lrho, mlrho2, s2.lcom, mcom2);
                if (mcom1 + mb != mcom2 + mlrho2 ||
                    abs(mcom1 + mb) > s2.lr ||
                    abs(mcom2 + mlrho2) > s2.lr) continue;

                res += P1 * cg(s2.lr, -mcom2-mlrho2, s2.lrho, mlrho2, s2.lcoup, -mcom2) *
                  cg(s2.lr, -mcom1-mb, b, mb, s1.lcoup, -mcom1) *
                  cg(s2.lcoup, -mcom2, s2.lcom, mcom2, ltotal, 0) *
                  cg(s1.lcoup, -mcom1, s1.lcom, mcom1, ltotal, 0);
              }
            }
          }
        }
        cache_G_rho.push_back(res);
      }
    }
  }
}

void angular_subsystem::generate_H(){

  cache_G_r.reserve(1);
  cache_G_rho.reserve(1);

  for (auto &s1 : states){
    for (auto &s2 : states_transformed){
      for (int kappa=0; kappa<lim_kappa; kappa++){
        double res = 0.0;
        for (int a=0; a<lim_special; a++){
          double E1 = E(kappa, kappa, 0, s1.lr, s1.lrho, s1.lcoup, a, s2.lrho, s1.lcoup);
          double G1 = G_r(s1, s2, a);
          res += E1 * G1;
        }
        cache_H_r.push_back(res);

        res = 0.0;
        for (int b=0; b<lim_special; b++){
          double E2 = E(kappa, kappa, 0, s1.lr, s1.lrho, s1.lcoup, s2.lr, b, s1.lcoup);
          double G2 = G_rho(s1, s2, b);
          res += E2 * G2;
        }
        cache_H_r.push_back(res);
      }
    }
  }
}

void angular_subsystem::print_caches(){
  print_array_to_file(cache_cg, name + ".cache_cg", "");
  print_array_to_file(cache_E, name + ".cache_E", "");
  print_array_to_file(cache_P, name + ".cache_P", "");
  print_array_to_file(cache_Q, name + ".cache_Q", "");
  print_array_to_file(cache_F, name + ".cache_F", "");
  print_array_to_file(cache_G_r, name + ".cache_G_r", "");
  print_array_to_file(cache_G_rho, name + ".cache_G_rho", "");
  print_array_to_file(cache_H_r, name + ".cache_H_r", "");
  print_array_to_file(cache_H_rho, name + ".cache_H_rho", "");
}
