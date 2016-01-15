#ifndef UTILITY_H
#define UTILITY_H

#include "general.h"
#include <time.h>

string extract_input(string input_string, int offset);
double extract_double(string var_name, string input_string, double default_value);
int extract_int(string var_name, string input_string, int default_value);
bool extract_bool(string var_name, string input_string, bool default_value);
string extract_string(string var_name, string input_string, string default_value);
bool extract_stringNew(string &output_string, string input_string, unsigned int &offset);
vector<double> extract_vector_double(string inputString, string varName, vector<double> default_value);
vector<int> extract_vector_int(string inputString, string varName, vector<int> default_value);

void print_matrix_to_file(double * &M, int len, string file_name);

template <typename T>
void print_vector(vector<T> V){
  for (unsigned int i = 0; i < V.size(); i++)
    cout << V[i] << endl;
  cout << endl;
}


void report_time();

void print_array_to_file(vector<double> &A, string file_name, string identifier);
bool read_array_from_file(vector<double> &tar, string file_name, string identifier, int desired_array_size);

#endif
