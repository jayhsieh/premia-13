#ifndef  _BS2D_STD2D_H
#define _BS2D_STD2D_H

#include "bs2d/bs2d.h"
#include "std2d/std2d.h"

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_random.h"
#include "numfunc.h"
#include "transopt.h"
#include "math/linsys.h"
#include <float.h>

/* 
static int Delta_Operator(double u1,double d1,double u2,double d2,double stock1,double stock2,double puu, double pud,double pdu, double pdd,double *ptdelta1,double *ptdelta2)
{
  *ptdelta1=((d2-1.)*(pdu-puu)+(u2-1.)*(pud-pdd))/(stock2*(u1-d1)*(u2-d2));
  *ptdelta2=((d1-1.)*(pud-puu)+(u1-1.)*(pdu-pdd))/(stock1*(u1-d1)*(u2-d2));
	
  return OK;
}
*/

extern int MOD_OPT(Delta_Operator)(double u1,double d1,double u2,double d2,double stock1,double stock2,double puu, double pud,double pdu, double pdd,double *ptdelta1,double *ptdelta2);



#endif
