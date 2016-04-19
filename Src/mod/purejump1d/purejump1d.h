#ifndef _PUREJUMP1D_H
#define _PUREJUMP1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD PUREJUMP1D

/*1D Pure Jump World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Mu;
  VAR Sigma;
  /*VAR Divid;*/
  VAR R;    
  VAR Beta;  
  VAR Nu;
} TYPEMOD;

#endif
