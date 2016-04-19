#ifndef _LMM1D_CGMY_H
#define _LMM1D_CGMY_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD LMM1D_CGMY

/*1D Libor Market Model CGMY Driven SDE*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR l0;
  VAR Sigma;
  VAR C;  
  VAR G;  
  VAR M;
  VAR Y;
  
} TYPEMOD;



#endif
