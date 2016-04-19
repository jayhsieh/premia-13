#ifndef _BS1D_H
#define _BS1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD BS1D

/*1D BlackScholes World*/     

typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Mu;
  VAR Sigma;
  VAR Divid;
  VAR R;       
} TYPEMOD;



#endif
