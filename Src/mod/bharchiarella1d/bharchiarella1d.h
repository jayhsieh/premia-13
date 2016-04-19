#ifndef _BharChiarella1D_H
#define _BharChiarella1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD BharChiarella1D

/*1D Bhar Chairella World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR alpha0;
  VAR alphar;
  VAR alphaf;
  VAR gamm;
  VAR tau;
  VAR lambda;
  VAR beta0;
  VAR beta1;
  VAR eta;   
} TYPEMOD;

#endif
