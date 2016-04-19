#ifndef _TIMEHES1D_H
#define _TIMEHES1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD TIMEHES1D

/*1D TIME DEPENDENT HESTON World*/    
 
typedef struct TYPEMOD{
  VAR T;
  VAR S0;
  VAR Divid;
  VAR R; 
  VAR Sigma0;
  VAR MeanReversion;
  VAR TimeDepParameters;
  VAR TimeStep;
} TYPEMOD;


#endif
