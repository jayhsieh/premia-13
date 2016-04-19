#ifndef _Affine3D_H
#define _Affine3D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD Affine3D

/*3D Affine World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR x01;
  VAR x02;
  VAR x03;
  VAR k1;
  VAR k2;
  VAR k3;
  VAR Sigma1;
  VAR Sigma2;
  VAR Sigma3;
  VAR shift;
  VAR Rho12;  
  VAR Rho13;  
  VAR Rho23; 
} TYPEMOD;


#endif
