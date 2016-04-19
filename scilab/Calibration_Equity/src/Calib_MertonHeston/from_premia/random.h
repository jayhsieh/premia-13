#ifndef __RANDOM_H__
#define __RANDOM_H__
#include "mathtools.h"

/*RandomGenerators*/
typedef struct Generator
{
  Label Name;
  void (*Compute)(int, double *);
  int RandOrQuasi;
  int Dimension;
} Generator;

int InitGenerator(int type_generator, int simulation_dim,long samples);
void Display_Generators(void);

int Rand_Or_Quasi(int type_generator);
int test_sobol(int type_generator);

double GaussQMC(int dimension, int create_or_retrieve, int index, int type_generator);
double GaussMC(int dimension, int create_or_retrieve, int index, int type_generator);

double premia_rand_normal(int generator);
double Gauss_BoxMuller(int generator);
int Bernoulli(double p, int generator);
long Poisson(double lambda, int type_generator);
double Uniform(int generator);
double D_Uniform(int dimension, int create_or_retrieve, int index, int type_generator, int mc_or_qmc);

double Inverse_erf(double x);
double Inverse_erf_Moro(double x);

extern double (*Gaussians[])(int,int,int,int);
int InitMc(void);

#endif
