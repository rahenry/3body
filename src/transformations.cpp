#include "three_body_system.h"
#include "states.h"
#include "gsl/gsl_blas.h"
#include "utility.h"

void three_body_system::calc_transformation(vector<double> &tar, state_full &s){
  bool quiet = 1;
  tar = vector<double>(9, 0.0);

  vector<double> bt;
  get_single_particle_transformation(bt, s.channel);

  gsl_matrix_view A = gsl_matrix_view_array(U_matrix_base->data, 3, 3);
  gsl_matrix_view B = gsl_matrix_view_array(&bt[0], 3, 3);
  gsl_matrix_view C = gsl_matrix_view_array(Uinv_matrix_base->data, 3, 3);
  gsl_matrix_view res = gsl_matrix_view_array(&tar[0], 3, 3);
  gsl_matrix *temp1 = gsl_matrix_alloc(3,3);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &A.matrix, &B.matrix, 0.0, temp1);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp1, &C.matrix, 0.0, &res.matrix);

  gsl_matrix_free(temp1); 

  if (!quiet){
    cout << "Transformation calculated: " << endl;
    print_vector(tar);
  }

}

void three_body_system::calc_transformation_coefficients(vector<double> &tar, vector<double> &t){

  bool quiet = 1;

  tar = vector<double>(6, 0.0);

  gsl_matrix *test_state_r = gsl_matrix_alloc(2,2);
  gsl_matrix_set_zero(test_state_r);
  gsl_matrix_set(test_state_r, 0, 0, 1.0);
  gsl_matrix *test_state_rho = gsl_matrix_alloc(2,2);
  gsl_matrix_set_zero(test_state_rho);
  gsl_matrix_set(test_state_rho, 1, 1, 1.0);

  gsl_matrix *transformation = gsl_matrix_alloc(2,2);
  gsl_matrix_set(transformation, 0, 0, t[0]);
  gsl_matrix_set(transformation, 0, 1, t[1]);
  gsl_matrix_set(transformation, 1, 0, t[3]);
  gsl_matrix_set(transformation, 1, 1, t[4]);

  gsl_matrix *temp1 = gsl_matrix_alloc(2,2);
  gsl_matrix *temp2 = gsl_matrix_alloc(2,2);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, test_state_r, transformation, 0.0, temp1);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, transformation, temp1, 0.0, temp2);

  tar[0] = gsl_matrix_get(temp2, 0, 0);
  tar[2] = gsl_matrix_get(temp2, 1, 0);
  tar[4] = gsl_matrix_get(temp2, 1, 1);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, test_state_rho, transformation, 0.0, temp1);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, transformation, temp1, 0.0, temp2);

  tar[1] = gsl_matrix_get(temp2, 0, 0);
  tar[3] = gsl_matrix_get(temp2, 1, 0);
  tar[5] = gsl_matrix_get(temp2, 1, 1);

  if (!quiet){
    cout << "Transformation coefficients calculated:" << endl;
    cout << tar[0] << " " << tar[1] << endl;
    cout << tar[2] << " " << tar[3] << endl;
    cout << tar[4] << " " << tar[5] << endl;

  }
  gsl_matrix_free(temp1); gsl_matrix_free(temp2); gsl_matrix_free(test_state_r); gsl_matrix_free(test_state_rho);
}

void three_body_system::calc_transformed_sums(vector<double> &tar, state_full &s1, state_full &s2){
  bool quiet = 1;
  vector<double> t1, t2;
  calc_transformation(t1, s1);
  calc_transformation(t2, s2);

  vector<double> tc1, tc2;
  calc_transformation_coefficients(tc1, t1);
  calc_transformation_coefficients(tc2, t2);

  apply_transformation(tar, s1, s2, tc1, tc2);

  if (!quiet){
    cout << "Calculated transformed sums: " << s1.channel << " -> " << s2.channel << endl;
  }
}


void three_body_system::get_single_particle_transformation(vector<double> &tar, int channel){
  double index = -1;
  if (channel > 0) index = channel + 2;
  if (channel < 0) index = -channel - 1;
  if (abs(channel) > 3) index = -1;
  if (index == -1){
    cout << "INVALID CHANNEL" << endl;
    return;
  }

  tar = single_particle_transformations[index];
}

