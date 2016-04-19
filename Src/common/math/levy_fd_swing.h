
#ifndef LEVYFDSWING_H
#define LEVYFDSWING_H

#include <iostream>
#include <fstream>
#include <vector>
#include "progonka.h"
#include "numerics.h"
#include "levy.h"

class GridSwing
{
  double Al;
  double Ar;
  double dx;
  int N;
    
 public :

  GridSwing(const double dAl, const double dAr, const int dN);
  inline double x(double i) const {return Al+i*dx;}
};

double init_cond_swing(const double x, const double S0,
                 const double K, const int product);

double bound_cond_swing(const double x, const double S0, const double K,const double rebate,
		  const double ttm, const double r,
		  const int product, const int product_type);

int price2_swing(int am,const Levy_measure & measure, int product,
		      int product_type, double r, double divid,double S0,
		      double K, double rebate,double Al, double Ar,
		      int Nspace, double T, int Ntime,int Nu,int Nd,double r_period, double & price0, double & delta0);

#endif

