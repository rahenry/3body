#ifndef RNG_H
#define RNG_H

#include <random>

static std::uniform_real_distribution<double> unif(0.0, 1.0);
static std::random_device rand_dev;
static std::mt19937_64 rand_engine(rand_dev());

double rng();

#endif
