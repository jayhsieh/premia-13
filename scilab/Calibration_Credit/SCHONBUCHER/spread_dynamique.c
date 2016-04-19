#include "spread_dynamique.h"

double spread_CDS(double t,double T,int M,int nb,double r,double R,double **trans){
double pl=0;
double dl=0;
double A;
double B;
int l,m,i,j;
double ***p;
p=malloc(M*sizeof(double**));

for(i=0;i<M;i++){
p[i]=malloc(M*sizeof(double*));
}

for(i=0;i<M;i++){
 for(j=0;j<M;j++){
 p[i][j]=malloc(nb*sizeof(double));
 }
}

 int n=0;
 p=proba(t,T,M,nb,trans);


/** We use the closed formula gin the Schonbucher paper**/

for(l=0;l<nb;l++){
  A=0;
  B=0;
  for(m=0;m<M;m++){
   A+=(1-R)*trans[m][l]*p[n][m][l];
   B+=(M-m*(1-R))*p[n][m][l];
  }
  dl+=(T-t)/(nb-1)*A*exp(-r*l*(T-t)/(nb-1)); /**default leg**/
  pl+=(T-t)/(nb-1)*B*exp(-r*l*(T-t)/(nb-1)); /**payment leg**/
}
for(i=0;i<M;i++){
  
 for(j=0;j<M;j++){
  free(p[i][j]);
 }
}
for(i=0;i<M;i++){
  free(p[i]);
  
}
free(p);



 return 10000*(dl/pl);

}


double spread_CDO(double t,double T,int M,int nb,double r,double a,double b,double R,double **trans){
double pl=0;
double dl=0;
double A;
double B;
int l,m,i,j;

double ***p;
p=malloc(M*sizeof(double**));
for(i=0;i<M;i++){
p[i]=malloc(M*sizeof(double*));

}
for(i=0;i<M;i++){
 for(j=0;j<M;j++){
 p[i][j]=malloc(nb*sizeof(double));
 }
}

 int n=0;
 p=proba(t,T,M,nb,trans); /** matrix of probability transition**/
 int NL=(a*M)*1./(1-R);  /** NL(number lower) is lower number which the tranche CDO is impacted**/ 
 int NU=(b*M)*1./(1-R); /** NU (number upper) is the upper number which the tranche CDo is impacted**/



for(l=0;l<nb;l++){
  A=0;
  B=0;
  if(a==0) NL=-1;
  for(m=NL+1;m<NU+1;m++){
   A+=(1-R)*trans[m][l]*p[n][m][l];
   B+=(1-R)*(NU-m)*p[n][m][l];
  }
 
  
  dl+=(T-t)/(nb-1)*A*exp(-r*l*(T-t)/(nb-1));
  pl+=(T-t)/(nb-1)*B*exp(-r*l*(T-t)/(nb-1));
}



for(i=0;i<M;i++){
  
 for(j=0;j<M;j++){
  free(p[i][j]);

 }
}

for(i=0;i<M;i++){
  free(p[i]);
  
  
}
free(p);

if(pl==0) return 0.0000001;

 if(b>0.03){
 
  return(10000*dl/pl);

/** return the spread in bps**/

 }

 else if(b<=0.03) return(100*(dl-0.05*pl)/((b-a)*M*(1-R)));

 /**return the upfront **/
}