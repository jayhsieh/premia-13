#include "tryschon.h"



/** Simulation of the probability transition between 0 and t**/

double  ***proba(double t,double T,int M,int nb,double **trans){
int i,j,d;
double s;
int k;
double ***p;
double h1;
double A;


p=malloc(M*sizeof(double**));
for(i=0;i<M;i++){
 p[i]=malloc(M*sizeof(double*));
 
}
for(i=0;i<M;i++){
 for(j=0;j<M;j++){
 p[i][j]=malloc(nb*sizeof(double));
 }
}

/** Initialisation of probability to jump to i for j in the  (t t+d)   **/

 for(i=0;i<M;i++){
   for(j=0;j<M;j++){
     for(k=0;k<nb;k++){
      if(j==i) p[i][j][k]=1;
      else     p[i][j][k]=0;
   }
  }
 }

/** first we simulate the diagonal of the matrix **/
 for(i=0;i<M;i++){
  for(d=1;d<nb;d++){
   s=0;
   A=0;
   h1=(T-t)/((nb-1));
  
  for(k=0;k<d;k++){
    A+=trans[i][k]*h1;
  } 
  p[i][i][d]=exp(-A);
 }
}
/** We use the recursion of Schonbucher to simulate the others coefficients of the matrix see the paper of Schonbucher**/
 for(i=0;i<M;i++){
  for(j=1;j<M;j++){
   
  if(j>i){
   for(d=1;d<nb;d++){
    h1=(T-t)/((nb-1)); 
    A=0;
    for(k=0;k<d;k++){
     A+=p[j][j][d]*(trans[j-1][k]*p[i][j-1][k]/p[j][j][k])*h1;
    }
    
    p[i][j][d]=A;
    
   
   }
  }
 }
}

 return p;

 for(i=0;i<M;i++){
  for(j=0;j<M;j++){
   free(p[i][j]);
  }
 }
 for(i=0;i<M;i++){
   free(p[i]);
 }
 free(p);

}



int simul_perte(double t,double T,int M,int nb,double **trans){

/** the lost L_t jump a time t with the intensity a_L*dt **/

int i,n=100;
int  L=0;;
double a=0.;


  a=unif();
  if(a<trans[L][nb-1]*t/n) L+=1;

return L;
}


