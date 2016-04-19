
/// \file base.h
/// \brief numerical constant
/// \author  M. Ciuca (MathFi, ENPC)
/// \note (C) Copyright Premia 8 - 2006, under Premia 8 Software license
//
//  Use, modification and distribution are subject to the 
//  Premia 8 Software license

#ifndef _BASE_H
#define _BASE_H

#include <iostream>
//#define NDEBUG
#include <cassert>
#include <cstdlib>

const double TOLERANCE  = 1.e-03;
const int    SIMPSON_NO = 20;
const double RIEMANN_NUM_INTEGR_PRECISION = 1.e-03;
const double INC = 1.e-06;

//inline double SQR(double x)
//{
// return x*x;
//}

//using namespace std;

struct DateRate
{
	//public:
	double date;
	double rate;
	DateRate(double _date=0, double _rate=1):
	date(_date), rate(_rate)
	{
		if((date < 0) || (rate <= 0))
		{
			std::cout << "DateRate Constructor: Incorrect input data!\n";
			exit(1);
		}

	}  
};





#endif

