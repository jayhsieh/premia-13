#ifndef _MER1D_H
#define _MER1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"
#include "enums.h"

#define TYPEMOD MER1D

/*1D Merton World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Mu;
  VAR Sigma;
  VAR Divid;
  VAR R;    
  VAR Lambda;  
  VAR Mean;
  VAR Variance;
} TYPEMOD;

#endif
