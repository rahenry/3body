#include "rng.h"

double rng(){
  return unif(rand_engine);
}

double random_double(double lower, double upper){
  return lower + rng() * (upper - lower);
}

int random_int(int lower, int upper){
  return lower + floor( rng() * (upper-lower+1) );


}

