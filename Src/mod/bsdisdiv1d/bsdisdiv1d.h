
#ifndef _BSDISDIV1D_H
#define _BSDISDIV1D_H
#include "optype.h"
#include "var.h"

#define TYPEMOD BSDISDIV1D

/*1D BlackScholes World with Discrete Dividends*/
typedef struct TYPEMOD{
  VAR Size;
  VAR T;
  VAR S0;
  VAR R;
  VAR Mu;
  VAR Sigma;
  VAR Amounts;
  VAR Dates;       
} TYPEMOD;


#endif
