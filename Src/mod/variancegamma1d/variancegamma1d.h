#ifndef _VARIANCEGAMMA1D_H
#define _VARIANCEGAMMA1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD VARIANCEGAMMA1D

/*1D VarianceGamma World*/     
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
