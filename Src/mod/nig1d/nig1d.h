#ifndef _NIG1D_H
#define _NIG1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"
#include "enums.h"

#define TYPEMOD NIG1D

/*1D NIG World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Mu;
  VAR Divid;
  VAR R;    
  VAR Sigma;  
  VAR Theta;
  VAR Kappa;
} TYPEMOD;


#endif
