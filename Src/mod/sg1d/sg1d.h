#ifndef _SquaredGaussian1D_H
#define _SquaredGaussian1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD SG1D

/*1D  SquaredGaussian World*/
typedef struct TYPEMOD{ 
  VAR T;
  VAR flat_flag;
  VAR a;
  VAR Sigma;
} TYPEMOD;

extern double MOD(GetYield)(TYPEMOD *pt);
extern char* MOD(GetCurve)(TYPEMOD *pt);

#endif
