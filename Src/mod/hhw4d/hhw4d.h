#ifndef _HHW4D_H
#define _HHW4D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD HHW4D

/*HHW4D World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR r;
  VAR MeanReversion;
  VAR x0;
  VAR kappa;
  VAR sigma; 
  VAR rho;
} TYPEMOD;

#endif
