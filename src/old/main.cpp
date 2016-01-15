#include "general.h"
#include "utility.h"
#include "three_body_system.h"

int main(int argc, char* argv[]){
  string input_file;
  GLOBAL_TIME = clock();

  int mode;
  if (argc == 1) input_file = "none";
  else input_file = string(argv[1]);

  three_body_system tbs(input_file);

  tbs.output_coordinate_matrices();

  angular_subsystem ang(&tbs);
  print_vector(ang.states);
  print_vector(ang.states_transformed);
  ang.print_caches();

  return 0;
}

  
