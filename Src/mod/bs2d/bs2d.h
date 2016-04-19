#ifndef  _BS2D_H
#define _BS2D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD BS2D

typedef struct TYPEMOD{ 
  VAR T;
  VAR S01;
  VAR Mu1;
  VAR Sigma1;
  VAR Divid1;
  VAR S02;
  VAR Mu2;
  VAR Sigma2;
  VAR Divid2;
  VAR Rho;
  VAR R;       
} TYPEMOD;


#endif
