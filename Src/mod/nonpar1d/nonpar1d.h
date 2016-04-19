#ifndef _NonPar1D_H
#define _NonPar1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD NONPAR1D

/*1D Non Parametric World*/     
typedef struct TYPEMOD{
  VAR S0;
  VAR R;   
  VAR implied_volatility;
  
} TYPEMOD;


#endif
