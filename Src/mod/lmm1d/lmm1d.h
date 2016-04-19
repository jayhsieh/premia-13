#ifndef _LMM1D_H
#define _LMM1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"
#include "enums.h"

#define TYPEMOD LMM1D

/*1D Libor Market Model World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR NbFactors;
  VAR l0;
  VAR Sigma;
  
} TYPEMOD;

#endif
