#ifndef _Vasicek1D_H
#define _Vasicek1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD Vasicek1D

/*1D Vasicek World*/     

typedef struct TYPEMOD{ 
  VAR T;
  VAR r0;
  VAR k;
  VAR Sigma;
  VAR theta;
} TYPEMOD;

#endif
