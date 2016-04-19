#ifndef _JarrowYildirim1D_H
#define _JarrowYildirim1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD JarrowYildirim1D

/*1D  Jarrow Yildirim World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR I0;
  VAR an;
  VAR ar;
  VAR sigman;
  VAR sigmar;
  VAR sigma_cpi;
  VAR Rhonr;
  VAR Rhoncpi;
  VAR Rhorcpi;
} TYPEMOD;

#endif
