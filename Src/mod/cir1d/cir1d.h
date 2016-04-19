#ifndef _Cir1D_H
#define _Cir1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD Cir1D

/*1D Cir World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR r0;
  VAR k;
  VAR Sigma;
  VAR theta;
} TYPEMOD;


#endif
