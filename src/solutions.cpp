#include "three_body_system.h"
#include "angular_subsystem.h"
#include "states.h"
#include "gsl/gsl_blas.h"
#include "controller.h"
#include "utility.h"
#include "eigenvalues.h"

void three_body_system::regenerate_matrices(){

  basis_size = states.size();
  S.resize(basis_size * basis_size);
  H.resize(basis_size * basis_size);
  vector<double> current_contributions;
  int current_index;

  for (auto & s1 : states){
    for (auto & s2: states){
      current_index = s1.index + basis_size * s2.index;
      calc_matrix_elements(current_contributions, s1, s2);
      S[current_index] = current_contributions[0];
      H[current_index] = 0.0;
      for (unsigned int i=1; i<current_contributions.size(); i++) H[current_index] += current_contributions[i];
    }
  }
}

void three_body_system::diagonalise(){

  double * H_temp = new double[basis_size*basis_size];
  double * S_temp = new double[basis_size*basis_size];

  //cout << "basis_size = " << basis_size << ", matrix_array_size = " << matrix_array_size << endl;
  for (int i=0; i<basis_size; i++){
    for (int j=0; j<basis_size; j++){
      //if (i > j) continue;
      int index1 = i + basis_size * j;
      S_temp[index1] = S[index1];
      H_temp[index1] = H[index1];
    }
  }

  evals.resize(basis_size);
  evecs.resize(basis_size * basis_size);

  //print_matrix_to_file(&S[0], basis_size, name+".S");
  //print_matrix_to_file(&H[0], basis_size, name+".H");
  eigs(&evals[0], &evecs[0], basis_size, true, H_temp, S_temp, basis_size);

  cout << "Energy eigenvalues: ";
  for (int i=0; i<cont->n_eigen; i++){
    cout << evals[i] << " ";
  }
  cout << endl;

  //evecs_test();

  delete [] S_temp;
  delete [] H_temp;

}
