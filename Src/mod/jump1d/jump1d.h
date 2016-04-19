#ifndef _JUMP1D_H
#define _JUMP1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD JUMP1D

/*1D Jump World for Swing Options*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Mu;
  VAR Sigma;
  VAR Divid;
  VAR R;    
  VAR Lambda;  
  VAR Mean;
} TYPEMOD;

#endif
