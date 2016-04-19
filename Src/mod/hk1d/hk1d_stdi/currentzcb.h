
#ifndef CURRENTZCB_H
#define CURRENTZCB_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double CurrentZCB (double T, int flat_flag, double r_flat, char* init);
/*computes P(0,T): if (flat_flag==0) then P(0,T)=exp(-r_flat*T) else
  P(0,T) is obtained from the datas in the file "init" and interpolation */



#endif     
