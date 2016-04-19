
/* cds_plini1.h + cds_plini1.cpp 7/10/05 
   CDS with deterministic Piecewise Linear Intensity 

  ::WARNING:: INTEREST RATES ARE FLAT (short rate is constant=r)
  


*/
#ifndef _INTENSITYCALIB_H
#define _INTENSITYCALIB_H

#include <string>

using namespace std;

#define PIECEWISE_NUMBER 5
#define SIMPSON_NO 20 // steps number in Simpson numerical integration



//******************************************************************************
//*******************         CDS CALIBRATION    *******************************
//******************************************************************************

int cds_plini_cali(double r, double Z,
				   int n, double *timesT,
				   int noCDS, 
				   double arrayCDS[][2],
				   double gamma[][2],
				   double f_tolerance = 0.0000001,
				   double tolerance = 0.0001,
				   int maxNoIterations = 50);
// given a spreads curve extract the (piecewise linear) deterministic
// default intensity				   

#endif
