#ifndef _DUPIRE1D_H
#define _DUPIRE1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD DUP1D

/*1D Dupire World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Mu;
  VAR Sigma;
  VAR Divid;
  VAR R;       
} TYPEMOD;


#endif
