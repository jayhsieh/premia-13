#ifndef _LMM_JUMP1D_H
#define _LMM_JUMP1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD LMM_JUMP1D

/*1D Libor Market Model World with Jumps*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR l0;
  VAR Sigma;
 
} TYPEMOD;

#endif
