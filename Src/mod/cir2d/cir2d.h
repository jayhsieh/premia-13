#ifndef _Cir2D_H
#define _Cir2D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD Cir2D

/*2D Cir World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR x01;
  VAR x02;
  VAR k1;
  VAR k2;
  VAR Sigma1;
  VAR Sigma2;
  VAR theta1;
  VAR theta2;
  VAR shift;
} TYPEMOD;


#endif
