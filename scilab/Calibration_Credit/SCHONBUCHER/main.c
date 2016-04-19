#include "stdio.h"
#include "spread_calib.h"

int main(){
int  i;
int M=125;           /** Number of firms**/
int nb=50;           /** Number of discrete points**/ 
double T=3;          /** Maturity of the spread**/
double t=1;          /**final time of the simulation of the path of the spread**/
double r=0.01;       /** interest rate**/
double *spread;      /** to get the spread for calibration**/
double vol=0.001;    /** Volatility of the transition rate**/
double *dym_sp;      /** to get the calibrate spread **/ 
double a=0.0;       /** (a,b) tranche of CDO**/
double b=0.03;
double R=0.4;        /** Recovery rate**/
int l1,l2,l3,m;

dym_sp=malloc((nb+1)*sizeof(double));

spread=malloc(24*sizeof(double));

/** Spread_CDS**/
spread[0]=13;        /** 3 Years**/
spread[1]=24;        /** 5 Years**/
spread[2]=33;        /** 7 Years**/
spread[3]=44;        /** 10 Years**/

/** Spread_CD0 [0,0.03]**/
spread[4]=-1.76;
spread[5]=11.25;
spread[6]=25.75;
spread[7]=40.5;
/**Spread_CDO [0.03,0.06]**/
spread[8]=0.4;
spread[9]=58;
spread[10]=135;
spread[11]=340;
/**Spread [0.06,0.09]**/
spread[12]=0.01;
spread[13]=0.15;
spread[14]=0.39;
spread[15]=0.99;
/**Spread [0.09,0.12]**/
spread[16]=0.000001;
spread[17]=0.07;
spread[18]=0.18;
spread[19]=0.46;
/**Spread [0.12,0.22]**/
spread[20]=0.000000001;
spread[21]=0.03;
spread[22]=0.06;
spread[23]=0.14;






/** The user can simulate the path of the CDO or the path of the index CDS **/

dym_sp=calib_spread_CDO(t,T,M,nb,vol,r,a,b,R,spread);

/**dym_sp=calib_spread_CDS(t,T,M,nb,vol,r,R,spread);**/


printf("***********************************************\n");

/**  print the path of the spread spread**/

printf("Path of the spread ");
printf("\n");

for(i=0;i<nb;i++){
printf("%f\n",dym_sp[i]);
} 



 return 0;

} 
