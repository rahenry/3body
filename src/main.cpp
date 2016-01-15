#include "general.h"
#include "utility.h"
#include "three_body_system.h"
#include "tests.h"
#include "controller.h"

int main(int argc, char* argv[]){
  string input_file;
  GLOBAL_TIME = clock();
  //rng_test();

  if (argc == 1) input_file = "none";
  else if (argc == 2) input_file = string(argv[1]);

  ifstream input_stream(input_file);
  stringstream buffer;
  buffer << input_stream.rdbuf();
  input_stream.close();
  int mode = int(extract_double(buffer.str(), "mode", 0));

  if (mode == -4){
    test4();
  }
  if (mode == -5){
    test5();
  }
  if (mode == -6){
    test6();
  }
  if (mode == 5){
    hyper_test(input_file);
  }

  if (mode == 0){
    cout << "mode = 0" << endl;
    cout << "input_file = " << input_file << endl;
    controller c(input_file);
    c.run_all();
  }


  //test1();
  //test2();
  //test3();
  //test4();
  //test5();
  //tbs.output_coordinate_matrices();

  //angular_subsystem ang(&tbs);
  //print_vector(ang.states);
  //print_vector(ang.states_transformed);
  //ang.print_caches();

  return 0;
}

  
