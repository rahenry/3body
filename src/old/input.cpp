#include "three_body_system.h"
#include "utility.h"

void three_body_system::get_file_data(){
  ifstream input_stream(name);
  stringstream buffer;
  buffer << input_stream.rdbuf();
  file_data = buffer.str();
}

int three_body_system::get_int(string var_name, int default_value){
  return extract_int(file_data, var_name, default_value);
}

double three_body_system::get_double(string var_name, double default_value){
  return extract_double(file_data, var_name, default_value);
}

bool three_body_system::get_bool(string var_name, bool default_value){
  return extract_bool(file_data, var_name, default_value);
}

string three_body_system::get_string(string var_name, string default_value){
  return extract_string(file_data, var_name, default_value);
}

vector<double> three_body_system::get_vector_double(string var_name, vector<double> default_value){
  return extract_vector_double(file_data, var_name, default_value);
}

vector<int> three_body_system::get_vector_int(string var_name, vector<int> default_value){
  return extract_vector_int(file_data, var_name, default_value);
}
