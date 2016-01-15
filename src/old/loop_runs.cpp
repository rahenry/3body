#include "three_body_system.h"

void three_body_system::loop_runs(){

  n_runs = mass_ratios.size() * omega_ratios.size() * max(r0s_unlike.size(), r0s_like.size());

  int run_index = 0;
  for (auto &mass_ratio : mass_ratios){
    for (auto &omega_ratio : omega_ratios){
    }
  }


}

