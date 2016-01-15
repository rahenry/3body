#include "three_body_system.h"
#include "controller.h"
#include "angular_subsystem.h"
#include "matrix_operations.h"
#include "gsl/gsl_blas.h"

double three_body_system::trapping_element(int i, int j, int channel){
  double res = 0.0;
  int cind = channel_index(channel);
  for (int k=0; k<3; k++){
    res += masses[k] * omegas[k] * omegas[k] *
      gsl_matrix_get(Uinv[cind], k, i) *
      gsl_matrix_get(Uinv[cind], k, j);
  }
  return res;
}

double three_body_system::kinetic_element(int i, int j, int channel){
  double res = 0.0;
  int cind = channel_index(channel);
  for (int k=0; k<3; k++){
    res += gsl_matrix_get(U[cind], i, k) *
      gsl_matrix_get(U[cind], j, k) / masses[k];
  }
  return res;
}

double three_body_system::U_element(int i, int j){
  double res = 0.0;

  if (j > i+1) res = 0.0;
  else if (j == i+1) res = -1.0;
  else{
    double msum = 0.0;
    for (int k=0; k <= i; k++) msum += masses[k];
    res = masses[j] / msum;
  }

  return res;
}

double three_body_system::Uinv_element(int i, int j){
  double res = 0.0;

  if (i > j+1) res = 0.0;
  else if (j == 2) res = 1.0;
  else{
    double msum = 0.0;
    for (int k=0; k <= j; k++) msum += masses[k];

    if (i == j+1) res = -msum / (msum + masses[j+1]);
    else res = masses[j+1] / (msum + masses[j+1]);
  }

  return res;
}

void three_body_system::generate_coordinate_matrices(){
  U_matrix_base = gsl_matrix_alloc(3, 3);
  Uinv_matrix_base = gsl_matrix_alloc(3, 3);
  for (int i=0; i<3; i++) for (int j=0; j<3; j++){
    gsl_matrix_set(U_matrix_base, i, j, U_element(i,j));
    gsl_matrix_set(Uinv_matrix_base, i, j, Uinv_element(i,j));
  }

  U.resize(6);
  Uinv.resize(6);
  for (auto channel : channels_allowed){
    int cind = channel_index(channel);
    U[cind] = gsl_matrix_alloc(3,3);
    Uinv[cind] = gsl_matrix_alloc(3,3);
    
    vector<double> permutation;
    get_permutation_matrix(permutation, channel);
    gsl_matrix_view A = gsl_matrix_view_array(U_matrix_base->data, 3, 3);
    gsl_matrix_view Ainv = gsl_matrix_view_array(Uinv_matrix_base->data, 3, 3);
    gsl_matrix_view P = gsl_matrix_view_array(&permutation[0], 3, 3);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &A.matrix, &P.matrix, 0.0, U[cind]);
    matrix_inverse(&permutation[0], 3);
    P = gsl_matrix_view_array(&permutation[0], 3, 3);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &P.matrix, &Ainv.matrix, 0.0, Uinv[cind]);

    for (int i=0; i<3; i++) for (int j=0; j<3; j++){
      trapping[cind][i][j] = trapping_element(i, j, channel);
      kinetic[cind][i][j] = kinetic_element(i, j, channel);
    }

    scale_factor_r[cind] = pow(trapping[cind][0][0] / kinetic[cind][0][0], -0.25);
    scale_factor_rho[cind] = pow(trapping[cind][1][1] / kinetic[cind][1][1], -0.25);
    scale_factor_com[cind] = pow(trapping[cind][2][2] / kinetic[cind][2][2], -0.25);
  }
}

void three_body_system::output_coordinate_matrices(){
  cout << "Printing coordinate matrices... " << endl << flush;
  ofstream out(name + ".coord");

  for (auto channel : channels_allowed){
    int cind = channel_index(channel);
    out << "channel = " << channel << endl;
    out << "U" << endl;
    for (int i=0; i<3; i++){
      for (int j=0; j<3; j++){
	out << setw(5) << gsl_matrix_get(U[cind], i, j) << " ";
      }
      out << endl << endl;
    }

    out << "Uinv" << endl;
    for (int i=0; i<3; i++){
      for (int j=0; j<3; j++){
	out << setw(5) << gsl_matrix_get(Uinv[cind], i, j) << " ";

      }
      out << endl << endl;
    }

    out << "trapping" << endl;
    for (int i=0; i<3; i++){
      for (int j=0; j<3; j++){
	out << setw(5) << trapping[cind][i][j] << " ";

      }
      out << endl << endl;
    }


    out << "_________" << endl;
  }

  out << "U_matrix_base" << endl;
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      out << setw(5) << gsl_matrix_get(U_matrix_base, i, j) << " ";
    }
    out << endl;
  }

  out << "Uinv_matrix_base" << endl;
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      out << setw(5) << gsl_matrix_get(Uinv_matrix_base, i, j) << " ";
    }
    out << endl;
  }
}

void three_body_system::report_length_scales(){
  ofstream o(name + ".length");
  o << "Single-particle oscillator lengths: " << endl;
  for (int i=0; i<3; i++){
    o << pow(masses[i]*omegas[i], -0.5) << endl;
  }

  o << endl << "Jacobian oscillator lengths: " << endl;
  int cind = channel_index(3);
  for (int i=0; i<3; i++){
    o << pow(trapping[cind][i][i] / kinetic[cind][i][i], -0.25) / sqrt(2.0)<< endl;
  }

  o << endl << "Interaction lengths: " << endl;
  o << "r0_like = " << r0_like << endl;
  o << "r0_unlike = " << r0_unlike << endl;
  



}


