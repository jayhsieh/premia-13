#ifndef _FPS2D_H
#define _FPS2D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD FPS2D

/*2D Fouque Papanicolau Sircar World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Divid;
  VAR R; 
  VAR InitialSlow;
  VAR InitialFast;
  VAR SigmaSlow;
  VAR SigmaFast;
  VAR MeanReversionSlow;
  VAR MeanReversionFast;
  VAR LongRunVarianceSlow;
  VAR LongRunVarianceFast;
  VAR Rho1;  
  VAR Rho2;  
  VAR Rho12;  
} TYPEMOD;


#endif
