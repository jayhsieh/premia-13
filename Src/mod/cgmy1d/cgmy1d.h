#ifndef _CGMY1D_H
#define _CGMY1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD CGMY1D

/*1D CGMY  World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Mu;
  VAR Divid;
  VAR R;    
  VAR C;  
  VAR G;  
  VAR M;
  VAR Y;
} TYPEMOD;

#endif
