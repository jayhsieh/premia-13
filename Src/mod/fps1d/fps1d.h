#ifndef _FPS1D_H
#define _FPS1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD FPS1D

/*1D Fouque Papanicolau Sircar World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Divid;
  VAR R; 
  VAR Sigma0;
  VAR MeanReversion;
  VAR LongRunVariance;
  VAR Rho;  
  VAR SigmaF;
} TYPEMOD;


#endif
