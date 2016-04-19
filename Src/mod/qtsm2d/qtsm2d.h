#ifndef _QTSM2D_H
#define _QTSM2D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD QTSM2D

/*QTSM2D World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR x;
  VAR d0;
  VAR d;
  VAR theta;
   VAR GammaV;
  VAR SigmaV;
  VAR KappaV;
 
} TYPEMOD;


#endif
