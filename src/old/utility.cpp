#include "utility.h"

string extract_input(string input_string, int offset)
{
  string result;
  char c;
  unsigned int i = offset;
  while(i < input_string.size()){
    c = input_string[i];
    i++;
    if (c == '\n') break;
    result += c;
  }
  return result;
}

double extract_double(string input_string, string var_name, double default_value)
{
  size_t found;
  string result;
  string s = var_name + "=";

  found = input_string.find(s);
  if (found != string::npos)
  {
    result = extract_input(input_string, s.size() + found);
    if (result == string()) return default_value;
    else return atof(result.c_str());
  }
  else
    return default_value;
}

int extract_int(string input_string, string var_name, int default_value){
  return int(extract_double(input_string, var_name, double(default_value)));
}

bool extract_bool(string input_string, string var_name, bool default_value){
  return bool(extract_double(input_string, var_name, double(default_value)));
}

string extract_string(string input_string, string var_name, string default_value){
  size_t found;
  string result;
  string s = var_name + "=";

  found = input_string.find(s);
  if (found != string::npos)
  {
    result = extract_input(input_string, s.size() + found);
    return result;
  }
  else
    return default_value;
}

bool extract_stringNew(string &output_string, string input_string, unsigned int &offset){
  char c;
  bool result = false;
  output_string = string();
  while(offset < input_string.size()){
    c = input_string[offset];
    offset++;
    if (c == ' '){
      result = true;
      break;
    }
    else if (c == ';'){
      break;
    }
    output_string += c;
  }
  return result;
}

vector<double> extract_vector_double(string input_string, string var_name, vector<double> default_value){
  size_t found;
  string extractedString;
  string s = var_name + "=";
  bool flag = true;
  unsigned int offset;
  vector<double> result;

  found = input_string.find(s);
  if (found != string::npos)
  {
    offset = s.size() + found;
    while (flag){
      flag = extract_stringNew(extractedString, input_string, offset);
      result.push_back(atof(extractedString.c_str()));
    }
    return result;
  }
  else return default_value;
}

vector<int> extract_vector_int(string input_string, string var_name, vector<int> default_value){
  size_t found;
  string extractedString;
  string s = var_name + "=";
  bool flag = true;
  unsigned int offset;
  vector<int> result;

  found = input_string.find(s);
  if (found != string::npos)
  {
    offset = s.size() + found;
    while (flag){
      flag = extract_stringNew(extractedString, input_string, offset);
      result.push_back(int(atof(extractedString.c_str())));
    }
    return result;
  }
  else return default_value;
}

void report_time(){
  cout << double( clock() - GLOBAL_TIME) / (double)CLOCKS_PER_SEC<< endl;
  GLOBAL_TIME = clock();
}
      
void print_matrix_to_file(double * &M, int len, string file_name){
  ofstream outFile;
  outFile.precision (numeric_limits<double>::digits10 + 1);
  cout << "Printing to " << file_name << endl;
  outFile.open(file_name);
  for (int i=0; i<len; i++){
    for (int j=0; j<len; j++){
      outFile << M[i+j*len] << " ";
    }
    outFile << endl;
  }
  outFile.close();
}

void print_array_to_file(vector<double> &A, string file_name, string identifier){
  ofstream outFile;
  outFile.open(file_name);
  if (identifier != ""){
    outFile << identifier << endl;
  }
  outFile << "size = " << A.size() << endl;
  outFile.flags(ios::scientific);
  outFile.precision (numeric_limits<double>::digits10 + 1);
  for (auto & element : A){
    if (element == 0.0) outFile << "0" << endl;
    else if (element == 1.0) outFile << "1" << endl;
    else outFile << element << endl;
  }
  outFile.close();
}

bool readArrayFromFile(vector<double> &tar, string file_name, string identifier, int desiredArraySize){
  tar = vector<double>();
  ifstream inFile;
  inFile.open(file_name);
  string firstLine;
  inFile >> firstLine;
  if (firstLine != identifier){
    return false;
  }

  int secondLine;
  inFile >> secondLine;
  if (secondLine < desiredArraySize){
    return false;
  }

  double a;
  while (inFile){
    inFile >> a;
    tar.push_back(a);
  }
  tar.pop_back();

  return true;
}
