#ifndef _TEMPEREDSTABLE1D_H
#define _TEMPEREDSTABLE1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD TEMPEREDSTABLE1D

/*1D TEMPEREDSTABLE World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Mu;
  VAR Divid;
  VAR R;    
  VAR AlphaPlus;  
  VAR AlphaMinus;  
  VAR LambdaPlus;
  VAR LambdaMinus;
  VAR CPlus;
  VAR CMinus;
} TYPEMOD;

#endif
