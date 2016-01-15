#include "three_body_system.h"
#include "controller.h"
#include "utility.h"
#include "findV0.h"

void three_body_system::generate_spectrum(){
  
  if (!cont->use_spectrum) return;

  cout << "Generating spectrum... " << endl;
  double ainv_grad = (cont->ainv_max - cont->ainv_min) / double(1.0 + cont->ainv_n_points);

  ofstream spectrum_output(name + ".spectrum");

  for (int i=0; i < cont->ainv_n_points; i++){
    double ainv_current = cont->ainv_min + i * ainv_grad;
    V0_like = calc_V0(r0_like, ainv_current);
    V0_unlike = calc_V0(r0_unlike, ainv_current);
    regenerate_matrices();
    diagonalise();
    for (int j=0; j<cont->n_eigen; j++){
      spectrum_output << evals[j] << " ";
    }
    spectrum_output << endl;
  }




}
