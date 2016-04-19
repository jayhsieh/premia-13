#include "copulas.h"
#include            "cdo.h"

#include            "maths.h"
/**Partie 2 méthodes d'inversions et parametres implicites****/

/**intensité implicite déduite des données du marché sur le CDS on utilise la méthode de Newton pour pouvoir calculer le zero
contrat CDS */

double       intens1(const prod *produit,
                     const double s)
{
double xlow= 0.; 
double F;
double F1;
double F2;
int MAX_ITERATIONS = 10;
double ACCURACY = 0.001;
double x=xlow;
int i;
for ( i=0;i<MAX_ITERATIONS;i++){ 
F=spread_CDS(produit,x)-s;
if (abso(F)<ACCURACY) return (x);
F1=spread_CDS(produit,x+ACCURACY)-s;
F2=spread_CDS(produit,x-ACCURACY)-s;
x = x - (2*ACCURACY*F)/(F1-F2);
}
return (0); 
}

double        intens2(const prod *produit,
                      const double s)
{

int MAX_ITERATIONS=10;
double xl,xh;
double x1=0;
double x2=0.5;
double f1=prix_contrat_CDS(produit,x1,s);
double f2=prix_contrat_CDS(produit,x2,s);
double dx,f,rts,df,dxold,temp;
double ACCURACY = 0.00001;
int j;
if(f1*f2>0) return (0);
if(f1==0) return (x1);
if(f2==0) return (x2);

if(f1<0.0){
	xl=x1;
	xh=x2;
}
else{
	xh=x1;
	xl=x2;
}
rts=0.5*(x1+x2);
dxold=abso(x2-x1);
dx=dxold;
f=prix_contrat_CDS(produit,rts,s);
df=(prix_contrat_CDS(produit,rts+ACCURACY,s)-prix_contrat_CDS(produit,rts-ACCURACY,s))/(2*ACCURACY);
for(j=0;j<MAX_ITERATIONS;j++){
	if((((rts-xh)*df-f)*((rts-xl)*df-f)>0)||((abso(2.0*f)>abso(dxold*df)))){
	 dxold=dx;
	 dx=0.5*(xh-xl);
	 rts=xl+dx;
	 if(xl==rts) return (rts);
	}
	else{
	dxold=dx;
	dx=f/df;
	temp=rts;
	rts-=dx;
	if(temp==rts) return (rts);
	}
	if(abso(dx)<ACCURACY) return (rts);
	f=prix_contrat_CDS(produit,rts,s);
    df=(prix_contrat_CDS(produit,rts+ACCURACY,s)-prix_contrat_CDS(produit,rts-ACCURACY,s))/(2*ACCURACY);
    if(f<0) xl=rts;
	else    xh=rts;
}
return (0);
}

/**Méthode de Newton + dichotomie  pour deteminer les tranches-correlations **/

/**On suppose qu'on a les données sur les 5 tranches + l'index on cherche à determiner la correlation implicite sur 
chaque tranche par inversion de la formule su spread **/


double     *tranche_correl(const prod *produit,
                           const double *s1,
                           const double s2)
{
/** s1 tableau des spreads des 5 tranches ,s2 spread de l'index CDS **/

double *result; 
prod*  (*prods);
int i,k,j;
int MAX_ITERATIONS=10;
double ACCURACY,x1,x2,f1,f2,xl,xh, dx,f,rts,df,dxold,temp;
ACCURACY=0.000001;
double lambda=intens2(produit,s2);

//double lambda=(s2/10000)*1./(1-produit->recov);
/** on définit les cinq produits CDO dont on cherche la correl !**/

   prods=malloc(5*sizeof(prod));
   for(i=0;i<5;i++){
   prods[i]=malloc(sizeof(prod));
  
   }
   prods[0]->att=0.0;
  

   for(i=1;i<5;i++){
   prods[i-1]->det=prods[i]->att=0.03*i;
  }
   prods[4]->det=0.22;

  
for(i=0;i<5;i++){
prods[i]->maturite=produit->maturite;
prods[i]->rate=produit->rate;
prods[i]->recov=produit->recov;
prods[i]->nb=produit->nb;
prods[i]->nominal=produit->nominal;
}



   result=malloc(5 *sizeof(double));
  
   for(k=0;k<5;k++){

   x1=0.0; //la correl est compris entre 0 et 1
   x2=0.99;
   f1=spread_CDO(prods[k],x1,lambda)-s1[k];
   f2=spread_CDO(prods[k],x2,lambda)-s1[k];

 //  if(f1*f2>=0){
   
 //  if(abso(f1)>abso(f2))           result[k]=x2*x2;
 //  else if(abso(f1)<=abso(f2))     result[k]=x1*x1;
//}

  //else {   

   if(f1<0.0){
   xl=x1;
   xh=x2;
  }
  else if(f1>=0.0){
  xh=x1;
  xl=x2;
  }

  rts=0.5*(x1+x2);
  dxold=abso(x2-x1);
  dx=dxold;
  f=spread_CDO(prods[k],rts,lambda)-s1[k];
  df=(spread_CDO(prods[k],rts+ACCURACY,lambda)-spread_CDO(prods[k],rts-ACCURACY,lambda))/(2*ACCURACY);
  j=0;
  do{  
	j=j+1;
	if((((rts-xh)*df-f)*((rts-xl)*df-f)>0)||((abso(2.0*f)>abso(dxold*df)))){
	 dxold=dx;
	 dx=0.5*(xh-xl);
	 rts=xl+dx;
	 if(xl==rts) break ;
	}
	else{
	dxold=dx;
	dx=f/df;
	temp=rts;
	rts-=dx;
	if(temp==rts) break;
	}
	if(abso(dx)<ACCURACY)  break;
	f=spread_CDO(prods[k],rts,lambda)-s1[k];
        df=(spread_CDO(prods[k],rts+ACCURACY,lambda)-spread_CDO(prods[k],rts-ACCURACY,lambda))/(2*ACCURACY);
        if(f<0) xl=rts;
        else    xh=rts;
   
     }while(j<MAX_ITERATIONS);
 
     result[k]=rts*rts;
   }

return (result);


}

/** Base correlation voir papier Hull and White perfect copula **/


double         *Base_correl(const prod* produit,
                            const double *s1,
                            const double s2)
{

/** s1 sera un tableau dont les valeurs représenteront les spread des 6 tranches de CDO
 s2 sera le spread du CDS**/

int y;
double *result;

double x1,x2,f1,f2,xl,xh, dx,f,rts,df,dxold,temp,c;
int MAX_ITERATIONS = 10;
double ACCURACY = 0.0001;
prod* (*prods);
prod* (*nvprods);

double F,F1,F2,L,rho;

int i,k,j;

double lambda=intens2(produit,s2);
//double lambda=(s2/10000)*1./(1-produit->recov) voir JP Morgan;

/** on définit les cinq produits CDO dont on cherche la correl !**/

   prods=malloc(5*sizeof(prod));
   for(i=0;i<5;i++){
   prods[i]=malloc(sizeof(prod));
   
 }

   prods[0]->att=0.0;
  

   for(i=1;i<5;i++){
   prods[i-1]->det=prods[i]->att=0.03*i;
  }
   prods[4]->det=0.22;

   

for(i=0;i<5;i++){
prods[i]->maturite=produit->maturite;
prods[i]->rate=produit->rate;
prods[i]->recov=produit->recov;
prods[i]->nb=produit->nb;
prods[i]->nominal=produit->nominal;
}

/** on doit définir quatres differents types de prod :[0,6],[0,9],[0,12],[0,22]***/

nvprods=malloc(5*sizeof(prod));
for(i=0;i<4;i++){
   nvprods[i]=malloc(sizeof(prod));
}

 for(i=0;i<4;i++){   
 nvprods[i]->att=0.0;
   
   if (i!=3)           nvprods[i]->det=0.03*(i+2);
   else if (i==3)      nvprods[i]->det=0.22;
 }


for(i=0;i<4;i++){
nvprods[i]->maturite=produit->maturite;
nvprods[i]->rate=produit->rate;
nvprods[i]->recov=produit->recov;
nvprods[i]->nb=produit->nb;
nvprods[i]->nominal=produit->nominal;
}


   result=malloc(5 *sizeof(double));
   result=tranche_correl(produit,s1,s2);
   
   for(i=0;i<5;i++){
   result[i]=sqrt(result[i]);
   }
   
   for(k=1;k<5;k++){
   x1=0.0; /**la correl est compris entre 0 et 1 **/
   x2=0.99;
   c=0;
   
   for(j=0;j<=k;j++){
   c+=loss(prods[j],result[j],lambda);
   
   }
   

   f1=loss(nvprods[k-1],x1,lambda)-c;
   f2=loss(nvprods[k-1],x2,lambda)-c;

//   if(f1*f2>=0){
   
  // if(abso(f1)>abso(f2))           result[k]=x2*x2;
   //else if(abso(f1)<=abso(f2))     result[k]=x1*x2;
   
   //}

  //else {   

   if(f1<0.0){
   xl=x1;
   xh=x2;
  }
  else if(f1>=0.0){
  xh=x1;
  xl=x2;
  }

  rts=0.5*(x1+x2);
  dxold=abso(x2-x1);
  dx=dxold;
  f=loss(nvprods[k-1],rts,lambda)-c;
  df=(loss(nvprods[k-1],rts+ACCURACY,lambda)-loss(nvprods[k-1],rts-ACCURACY,lambda))/(2*ACCURACY);
  j=0;
  do{  
	j=j+1;
	if((((rts-xh)*df-f)*((rts-xl)*df-f)>0)||((abso(2.0*f)>abso(dxold*df)))){
	 dxold=dx;
	 dx=0.5*(xh-xl);
	 rts=xl+dx;
	 if(xl==rts) break ;
	}
	else{
	dxold=dx;
	dx=f/df;
	temp=rts;
	rts-=dx;
	if(temp==rts) break;
	}
	if(abso(dx)<ACCURACY)  break;
	f=loss(nvprods[k-1],rts,lambda)-c;
        df=(loss(nvprods[k-1],rts+ACCURACY,lambda)-loss(nvprods[k-1],rts-ACCURACY,lambda))/(2*ACCURACY);
        if(f<0) xl=rts;
        else    xh=rts;
   
     }while(j<MAX_ITERATIONS);
 
     result[k]=rts*rts; 
    }
      

  result[0]=result[0]*result[0];

  return (result);

} 






















