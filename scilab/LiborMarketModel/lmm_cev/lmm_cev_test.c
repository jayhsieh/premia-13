#include "lmm_cev_pricer.h"
#include <stdio.h>
int main(int argc,char* argv[])
{
  double* out=malloc(sizeof(double)*2);

  cev_price(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),0,0.4, 0.06, 0.065,4.,10.,atoi(argv[4]),out);
  printf("%f\t%f\n",out[0],out[1]);
  free(out);
  return 0;  
}
