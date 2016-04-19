#include "spread_calib.h"


double *calib_spread_CDS(double t,double T,int M,int nb,double vol,double r,double R,double *spread){
double pl=0;
double dl=0;
double A;
double B;
int l,m,i,j;
double ***trans;
double *dym_spread;
int n;
dym_spread=malloc((nb+1)*sizeof(double));
double ***p;
double **Ar;
int d;

Ar=malloc(M*sizeof(double*));
p=malloc(M*sizeof(double**));

for(i=0;i<M;i++){
p[i]=malloc(M*sizeof(double*));
Ar[i]=malloc(nb*sizeof(double));
}

for(i=0;i<M;i++){
 for(j=0;j<M;j++){
 p[i][j]=malloc(nb*sizeof(double));
 }
}
trans=malloc(M*sizeof(double**));

 for(i=0;i<M;i++){
   trans[i]=malloc(nb*sizeof(double*));
 }
 for(i=0;i<M;i++){
   for(j=0;j<nb;j++){
     trans[i][j]=malloc((nb+1)*sizeof(double));
   }
  }
 trans=calib_rates_CDS(t,T,nb,vol,r,spread);/** we get  the path of the transition rate**/  
 


for(d=0;d<nb;d++){
 n=trans[0][d][nb];
 
 for(i=0;i<M;i++){
  for(j=0;j<nb;j++){
   Ar[i][j]=trans[i][j][d];
  }
 }

 p=proba(d*t/(nb-1),T,M,nb,Ar);/** we get the probability transition for each step**/
 dl=0;
 pl=0;
 for(l=0;l<nb;l++){
  A=0;
  B=0;
  for(m=0;m<M;m++){
   A+=(1-R)*Ar[m][l]*p[n][m][l];
   B+=(M-m*(1-R))*p[n][m][l]; 
  }
  dl+=(T-t*d/(nb-1))/(nb-1)*A*exp(-r*l*(T-t*d/(nb-1))/(nb-1)); /** payment leg**/
  pl+=(T-t*d/(nb-1))/(nb-1)*B*exp(-r*l*(T-t*d/(nb-1))/(nb-1)); /**default leg**/
  if(pl==0) pl=0.00000001;
  dym_spread[d]=10000*dl/pl;
 }
}

for(i=0;i<M;i++){
  for(j=0;j<nb;j++){
    free(trans[i][j]);
    
  }
}

for(i=0;i<M;i++){
  
 for(j=0;j<M;j++){
  free(p[i][j]);
  
 }
}
for(i=0;i<M;i++){
  free(p[i]);
  free(trans[i]);
  free(Ar[i]);
}
free(p);
free(trans);
free(Ar);

 dym_spread[nb]=n;

 return dym_spread;

 free(dym_spread);
}
 
double *calib_spread_CDO(double t,double T,int M,int nb,double vol,double r,double a,double b,double R,double *spread){
double pl=0;
double dl=0;
double A;
double B;
int l,m,i,j;
double ***trans;
double *dym_spread;
int n;
dym_spread=malloc((nb+1)*sizeof(double));
double ***p;
double **Ar;
int d;

Ar=malloc(M*sizeof(double*));
p=malloc(M*sizeof(double**));

for(i=0;i<M;i++){
p[i]=malloc(M*sizeof(double*));
Ar[i]=malloc(nb*sizeof(double));
}

for(i=0;i<M;i++){
 for(j=0;j<M;j++){
 p[i][j]=malloc(nb*sizeof(double));
 }
}
trans=malloc(M*sizeof(double**));

 for(i=0;i<M;i++){
   trans[i]=malloc(nb*sizeof(double*));
 }
 for(i=0;i<M;i++){
   for(j=0;j<nb;j++){
     trans[i][j]=malloc((nb+1)*sizeof(double));
   }
  }

trans=calib_rates_CDO(t,T,nb,vol,r,spread);

int NL=(a*M)*1./(1-R);
int NU=(b*M)*1./(1-R);

 
for(d=0;d<nb;d++){
 n=trans[0][0][nb];
 
 for(i=0;i<M;i++){
  for(j=0;j<nb;j++){
   Ar[i][j]=trans[i][j][d];
  }
 }

 p=proba(d*t/(nb-1),T,M,nb,Ar);
 dl=0;
 pl=0;

 if(a==0) NL=-1;

 for(l=0;l<nb;l++){
  A=0;
  B=0;
  for(m=NL+1;m<NU+1;m++){
   A+=(1-R)*Ar[m][l]*p[n][m][l];
   B+=(1-R)*(NU-m)*p[n][m][l];
  }
  
  dl+=(T-t*d/(nb-1))/(nb-1)*A*exp(-r*l*(T-t*d/(nb-1))/(nb-1));
  pl+=(T-t*d/(nb-1))/(nb-1)*B*exp(-r*l*(T-t*d/(nb-1))/(nb-1));
 }
  if(b>0.03){
  if(pl==0) pl=0.00000001;
  dym_spread[d]=10000*dl/pl;
 }
 else if(b<=0.03) dym_spread[d]=100*(dl-0.05*pl)/((b-a)*M*(1-R));
  
/**upfront**/ 
 
}
for(i=0;i<M;i++){
  for(j=0;j<nb;j++){
    free(trans[i][j]);
   
  }
}

for(i=0;i<M;i++){
  
 for(j=0;j<M;j++){
  free(p[i][j]);
  
 }
}
for(i=0;i<M;i++){
  free(p[i]);
  free(trans[i]);
  free(Ar[i]);
}
free(p);
free(trans);
free(Ar);
dym_spread[nb]=n;

 return dym_spread;

free(dym_spread);
}

