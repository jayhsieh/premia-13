#ifndef _CirPlus2D_H
#define _CirPlus2D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD CirPP2D

/*2D Cir++ World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR flat_flag;
  VAR InitialYieldsR;
  VAR aR;
  VAR bR;
  VAR SigmaR;
  VAR InitialYieldsI;
  VAR aI;
  VAR bI;
  VAR SigmaI;
  VAR Rho;
} TYPEMOD;


#endif
