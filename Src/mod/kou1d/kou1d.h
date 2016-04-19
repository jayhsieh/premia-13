#ifndef _KOU1D_H
#define _KOU1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD KOU1D

/*1D Kou World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Mu;
  VAR Divid;
  VAR R;    
  VAR Sigma;  
  VAR Lambda;
  VAR LambdaPlus;
  VAR LambdaMinus;
  VAR P;
} TYPEMOD;

#endif
