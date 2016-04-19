#ifndef  _BS1D_STD_H
#define _BS1D_STD_H

#include "bs1d.h"
#include "std.h"

#include "mathtools.h"
#include "random.h"
#include "numfunc.h"
#include "transopt.h"
#include "linsys.h"

static double Nd1(double s,double r,double divid,double sigma,double T,double K) 
	{
	double d1=(log(s/K)+(r-divid+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
	
	return  N(d1); 
	}

#endif
