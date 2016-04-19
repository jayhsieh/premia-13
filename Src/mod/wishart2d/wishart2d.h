#ifndef _WISHART2D_H
#define _WISHART2D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD WISHART2D

/*WISHART2D World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR R;
  VAR Divid;
  VAR alpha;
  VAR b;
  VAR Sigma0; 
  VAR Q;
   
} TYPEMOD;


#endif
