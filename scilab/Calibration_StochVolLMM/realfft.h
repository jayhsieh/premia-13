#ifndef REALFFT
#define REALFFT

#include <vector>
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <complex>
using namespace std;
using std::string;
#define PI 3.142
const complex<double> I(0,1);


void realfastfouriertransform(double *a, int tnn, int inversefft);

#endif

