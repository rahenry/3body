#include "three_body_system.h"
#include "states.h"
#include "gsl/gsl_blas.h"
#include "utility.h"
#include "matrix_operations.h"

void three_body_system::get_permutation_matrix(vector<double> &tar, int channel){
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

void three_body_system::calc_coordinate_transformation(vector<double> &tar, int channel_initial, int channel_final){
  tar = vector<double>(9, 0.0);

  // Get the single-particle permutation matrices
  /*vector<double> permutation_initial, permutation_final;
  get_permutation_matrix(permutation_initial, channel_initial);
  get_permutation_matrix(permutation_final, channel_final);
  matrix_inverse(&permutation_initial[0], 3);

  // The result is U.Pf.Pi^-1.U^-1
  gsl_matrix_view U = gsl_matrix_view_array(U_matrix_base->data, 3, 3);
  gsl_matrix_view Pf = gsl_matrix_view_array(&permutation_final[0], 3, 3);
  gsl_matrix_view Pi = gsl_matrix_view_array(&permutation_initial[0], 3, 3);
  gsl_matrix_view Uinv = gsl_matrix_view_array(Uinv_matrix_base->data, 3, 3);

  gsl_matrix_view res = gsl_matrix_view_array(&tar[0], 3, 3);
  gsl_matrix *temp1 = gsl_matrix_alloc(3, 3);
  gsl_matrix *temp2 = gsl_matrix_alloc(3, 3);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &U.matrix, &Pf.matrix, 0.0, temp1);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Pi.matrix, &Uinv.matrix, 0.0, temp2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp1, temp2, 0.0, &res.matrix);

  gsl_matrix_free(temp1);
  gsl_matrix_free(temp2);*/

  int cind_initial = channel_index(channel_initial);
  int cind_final = channel_index(channel_final);
  gsl_matrix_view Uinv_old = gsl_matrix_view_array(Uinv[cind_initial]->data, 3, 3);
  gsl_matrix_view U_new = gsl_matrix_view_array(U[cind_final]->data, 3, 3);
  gsl_matrix_view res = gsl_matrix_view_array(&tar[0], 3, 3);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &U_new.matrix, &Uinv_old.matrix, 0.0, &res.matrix);

}

void three_body_system::calc_basis_transformation(vector<double> &tar, vector<double> &t){

  // calculates the transformation coefficients of the gaussian widths in a basis element
  // e.g. for applying a transformation X, the new alpha_r is X_transverse.[1 0; 0 0].X

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
  gsl_matrix_free(temp1); gsl_matrix_free(temp2); gsl_matrix_free(test_state_r); gsl_matrix_free(test_state_rho); gsl_matrix_free(transformation);
}

void three_body_system::apply_transformation(vector<double> &tar, state_full &s1, state_full &s2, vector<double> &tc1, vector<double> &tc2){
  // applies a transformation to the gaussian widths of two basis elements
  // tc1 and tc2 specify the transformation coefficients for the two basis elements
  // the result (which is three gaussian widths representing the transformed state) is stored in tar
  tar = vector<double>(3, 0.0);

  tar[0] += tc1[0] * s1.alpha_r;
  tar[0] += tc1[1] * s1.alpha_rho;
  tar[0] += tc2[0] * s2.alpha_r;
  tar[0] += tc2[1] * s2.alpha_rho;

  tar[1] += tc1[2] * s1.alpha_r;
  tar[1] += tc1[3] * s1.alpha_rho;
  tar[1] += tc2[2] * s2.alpha_r;
  tar[1] += tc2[3] * s2.alpha_rho;

  tar[2] += tc1[4] * s1.alpha_r;
  tar[2] += tc1[5] * s1.alpha_rho;
  tar[2] += tc2[4] * s2.alpha_r;
  tar[2] += tc2[5] * s2.alpha_rho;
}


void three_body_system::calc_transformed_state(vector<double> &tar, state_full &s1, state_full &s2, int channel_final){

  // takes two states and transforms to the specified channel
  // returns the gaussian widths of the final state (r, rho and r.rho)
  bool quiet = 1;
  vector<double> t1, t2;
  calc_coordinate_transformation(t1, s1.channel, channel_final);
  calc_coordinate_transformation(t2, s2.channel, channel_final);
  calc_coordinate_transformation(t1, channel_final, s1.channel);
  calc_coordinate_transformation(t2, channel_final, s2.channel);

  vector<double> tc1, tc2;
  calc_basis_transformation(tc1, t1);
  calc_basis_transformation(tc2, t2);

  apply_transformation(tar, s1, s2, tc1, tc2);

  if (!quiet){
    cout << "Calculated transformed sums: " << s1.channel << " -> " << s2.channel << endl;
  }
}

void three_body_system::initialise_transformations(){
  single_particle_transformations = vector< vector<double> >();

  single_particle_transformations.push_back(vector<double>{0, 0, 1, 0, 1, 0, 1, 0, 0});
  single_particle_transformations.push_back(vector<double>{1, 0, 0, 0, 0, 1, 0, 1, 0});
  single_particle_transformations.push_back(vector<double>{0, 1, 0, 1, 0, 0, 0, 0, 1});
  single_particle_transformations.push_back(vector<double>{0, 1, 0, 0, 0, 1, 1, 0, 0});
  single_particle_transformations.push_back(vector<double>{0, 0, 1, 1, 0, 0, 0, 1, 0});
  single_particle_transformations.push_back(vector<double>{1, 0, 0, 0, 1, 0, 0, 0, 1});
}


void three_body_system::generate_B_coef_caches(){
  cout << "Generating B caches... " << endl;
  vector<double> channel_start_allowed{-1, -2, -3, 1, 2, 3};
  vector<double> channel_end_allowed{-1, -2, -3, 1, 2, 3};

  for (auto channel_start : channel_start_allowed){
    B_coef_caches.push_back( vector<B_coef_cache>() );
    for (auto channel_end: channel_end_allowed){
      B_coef_caches.back().push_back(B_coef_cache(this, channel_start, channel_end));
    }
  }
}

int three_body_system::channel_index(int channel){
  int index;
  if (abs(channel) > 3 || channel==0) cout << "BAD CHANNEL" << endl;
  if (channel < 0) index= -channel - 1;
  else index = channel + 2;
  return index;
}

double three_body_system::B(state_full &state_start, int transformation_index, int channel_end){
  bool quiet = 1;
  int index_start, index_end;

  if (state_start.channel < 0) index_start = -state_start.channel - 1;
  else index_start = state_start.channel + 2;

  if (channel_end < 0) index_end = -channel_end - 1;
  else index_end = channel_end + 2;

  if (!quiet) cout << "Retrieving B coef (" << state_start.channel << ", " << channel_end << " | " << state_start.ang->index << ", " << transformation_index << " | " << index_start << " " << index_end << ")... " << flush;

  double res = B_coef_caches[index_start][index_end].B_coefs[state_start.ang->index][transformation_index];
  if (!quiet) cout << res << endl << flush;

  return res;
}
