#include "general.h"
#include "three_body_system.h"

double three_body_system::trapping_element(int i, int j){
  double res = 0.0;
  for (int k=0; k<3; k++){
    res += masses[k] * omegas[k] * omegas[k] * Uinv[k][i] * Uinv[k][j];
  }
  return res;
}

double three_body_system::kinetic_element(int i, int j){
  double res = 0.0;
  for (int k=0; k<3; k++){
    res += U[i][k] * U[j][k] / masses[k];
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
    for (int k=0; k <= i; k++) msum += masses[k];

    if (j == i+1) res = -msum / (msum + masses[j+1]);
    else res = masses[j+1] / (msum + masses[j+1]);
  }

  return res;
}

void three_body_system::generate_coordinate_matrices(){
  for (int i=0; i<3; i++) for (int j=0; j<3; j++){
    U[i][j] = U_element(i, j);
    Uinv[i][j] = Uinv_element(i, j);
  }

  for (int i=0; i<3; i++) for (int j=0; j<3; j++){
    trapping[i][j] = trapping_element(i, j);
    kinetic[i][j] = kinetic_element(i, j);
  }
}

void three_body_system::output_coordinate_matrices(){
  ofstream out(name + ".coord");

  out << "U" << endl;
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      printf( "%8.8f " , U[i][j]);
    }
    printf( "\n" );
  }
}



