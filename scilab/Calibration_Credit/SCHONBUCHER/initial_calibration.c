#include "initial_calibration.h"

double *initial_rates_CDS(double T,int n,double r,double *spread){
int i,j;
int l1,l2,l3,k;
double f,f1,sp,s1,s,s2;
int Max_iter=20;/** Maximum of the iteration in the Newton method**/ 
double eps=0.0000001;
double *init;
double **a;
int nb; 
 
int M=125; /** Number of firms for the itraxx this number is fixed**/
double R=0.4;/** We fix the recovery for the calibration**/

l1=(0.3*n);
l2=(0.5*n);
l3=0.7*n;

/** l1,l2,l3 is used to calibrate  a_n(0,s)  0<s<3 a_n(0,s)=cste1, 3<s<5 a_n(0,s)=cte2, 5<s<7 a_n(0,s)=cste3, 7<s<10 a_n(0,s)=cste4**/

 /** the size of init is a function of the maturity T**/

 if(T<=3){
   nb=1;
   
  }
  else{
     if(T<=5){
        nb=2;
        
      }
     else{
          if(T<=7){
             nb=3;
             
          }
          else{
           nb=4;
           
          }
     }
    }

init=malloc(nb*sizeof(double*));

a=malloc(M*sizeof(double));

for(i=0;i<M;i++){
 a[i]=malloc(n*sizeof(double));
}
s=0.1;

for(j=0;j<M;j++){
  for(i=0;i<n;i++){
   a[j][i]=s;
  }
}
/**Méthode de Newton pour calibrer les transition rates**/

sp=spread_CDS(0,3,M,n,r,R,a);
 
j=0;

/** We use Newton method to fit the spread with the good initial rate a(0,s)**/ 

do{
 f=sp-spread[0];
 j++;

for(k=0;k<M;k++){
 for(i=0;i<n;i++){
  a[k][i]=s+eps;
 }
}

 sp=spread_CDS(0,3,M,n,r,R,a);
 f1=sp-spread[0];
 s1=s-(f/(f1-f))*eps;
 if(fabs(s1-s)<eps) break;
 else s=s1;

}while(j<Max_iter);

init[0]=s1;

nb--;


if(nb>0){
 s=0.1;
 for(k=0;k<M;k++){
  for(i=0;i<l1;i++){
   a[k][i]=init[0];
  }
  for(i=l1;i<n;i++){
   a[k][i]=s;
  }
}

 sp=spread_CDS(0,5,M,n,r,R,a);
 j=0;

do{
  f=sp-spread[1];
  j++;
 for(k=0;k<M;k++){
  for(i=l1;i<n;i++){
   a[k][i]=s+eps;
  }
 }
 sp=spread_CDS(0,5,M,n,r,R,a);
 f1=sp-spread[1];
 s2=s-(f/(f1-f))*eps;
 if(fabs(s2-s)<eps) break;
 else s=s2;

 }while(j<Max_iter);
}
 init[1]=s2;
 
 nb--;

if(nb>0){
 s=0.1;
 
 for(k=0;k<M;k++){
  for(i=0;i<l1;i++){
   a[k][i]=init[0];
  }
  for(i=l1;i<l2;i++){
   a[k][i]=init[1];
  }
  for(i=l2;i<n;i++){
   a[k][i]=s;
  }
}
 

 sp=spread_CDS(0,7,M,n,r,R,a);
 j=0;

do{
  f=sp-spread[2];
  j++;
 
for(k=0;k<M;k++){ 
 for(i=l2;i<n;i++){
   a[k][i]=s+eps;
  }
 }
 sp=spread_CDS(0,7,M,n,r,R,a);
 f1=sp-spread[2];
 s2=s-(f/(f1-f))*eps;
 if(fabs(s2-s)<eps) break;
 else s=s2;

 }while(j<Max_iter);
}
 init[2]=s2;
 
 nb--;

if(nb>0){
 s=0.1;
 for(k=0;k<M;k++){
  for(i=0;i<l1;i++){
   a[k][i]=init[0];
  }
  for(i=l1;i<l2;i++){
   a[k][i]=init[1];
  }
  for(i=l2;i<l3;i++){
   a[k][i]=init[2];
  }
  for(i=l3;i<n;i++){
   a[k][i]=s;
  }
 }
 sp=spread_CDS(0,10,M,n,r,R,a);
 j=0;

do{
  f=sp-spread[3];
  j++;

for(k=0;k<M;k++){ 
 for(i=l3;i<n;i++){
   a[k][i]=s+eps;
  }
 }
 sp=spread_CDS(0,10,M,n,r,R,a);
 f1=sp-spread[3];
 s2=s-(f/(f1-f))*eps;
 if(fabs(f1)<eps) break;
 else s=s2;

 }while(j<Max_iter);
}
 init[3]=s2;


return init;

for(i=0;i<M;i++){
 free(a[i]);
}
free(a);

free(init);



}




/** this function return a_n(0,s); 0<s<T; 0<n<125 using the itraxx spread of CDO**/

double *initial_rates_CDO(double T,int n,double r,double *spread){

/**Spread est un vecteur composé de 6 vecteurs superposés Spread 1 représente l'évolution du spread de 3 à 10 ans, Spreadi i=2..6 représente l'évolution de la tranchei de CDO de 3 à 10 ans **/
/**On cherche dans cette étude à trouver les transition_rates initiaux tels qu'on approchera au mieux ces spreads par les formules fermées de Schonbucher**/

int i,j;
int l1,l2,l3,k;
double sp,s1,s,s2;
double f,f1;
int Max_iter=20;
double eps=0.0000001;
double *init;
double **a;


int m;
int nb; 
int M=125; /**The number of firms is fixed for the calibration**/
double R=0.4; /** the recovery is fixed for the calibration**/
 
l1=(0.3*n);
l2=(0.5*n);
l3=0.7*n;

 /** The size of the spread=5*nb is a function of the maturity T **/


 if(T<=3){
   nb=1;
   m=l1;
  }
  else{
     if(T<=5){
        nb=2;
        m=l2;
      }
     else{
          if(T<=7){
             nb=3;
             m=l3;
          }
          else{
           nb=4;
           m=n;
          }
     }
    }

a=malloc(M*sizeof(double*));

/** to return the initial rate**/

init=malloc((nb*5)*sizeof(double));


for(i=0;i<M;i++){
 a[i]=malloc(n*sizeof(double));
}

s=0.1;
for(j=0;j<M;j++){
 for(i=0;i<n;i++){
   a[j][i]=s;
  }
}

/**We use Newton method to fit the spread**/

/** first we find the good transition rate  using the first itraxx tranche (0,0.03) and ...**/

sp=spread_CDO(0,3, M,n,r,0.0,0.03,R,a); 

j=0;

do{
 f=sp-spread[4];
 j++;

for(k=0;k<7;k++){
 for(i=0;i<n;i++){
  a[k][i]=s+eps;
 }
}
sp=spread_CDO(0,3,M,n,r,0.0,0.03,R,a);
f1=sp-spread[4];
s1=s-(f/(f1-f))*eps;

 if(fabs(s1-s)<eps) break;
 else s=s1;

}while(j<Max_iter);

 init[0]=s1;


 j=0;

 s=0.01;
 for(j=0;j<M;j++){
  for(i=0;i<n;i++){
   a[j][i]=s;
  }
}
 
sp=(spread_CDO(0,3,M,n,r,0.03,0.06,R,a));
 
do{ 
   j++;
   f=sp-spread[8];
   for(k=7;k<13;k++){
    for(i=0;i<n;i++){
     a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,3,M,n,r,0.03,0.06,R,a));
 f1=sp-spread[8];
 s2=s-(f/(f1-f))*eps;
 
  if(fabs(s2-s)<eps) break;
 
  else s=s2;

 }while(j<Max_iter);
 
 init[1]=s;

 s=0.01;

 for(j=0;j<M;j++){
  for(i=0;i<n;i++){
   a[j][i]=s;
  }
 }
 j=0;
 
 sp=(spread_CDO(0,3, M,n,r,0.06,0.09,R,a));

do{ 
  
  j++;
  f=sp-spread[12];
  for(k=13;k<19;k++){
   for(i=0;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,3,M,n,r,0.06,0.09,R,a));
 f1=sp-spread[12];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

 }while(j<Max_iter);

  init[2]=s;

 s=0.01;

 for(j=0;j<M;j++){
  for(i=0;i<n;i++){
   a[j][i]=s;
  }
 }
 j=0;
 
 sp=(spread_CDO(0,3, M,n,r,0.09,0.12,R,a));
 
do{ 
   j++;
  f=sp-spread[16];
  for(k=19;k<26;k++){
   for(i=0;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,3,M,n,r,0.09,0.12,R,a));
 f1=sp-spread[16];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

 }while(j<Max_iter);

 init[3]=s;
 
 s=0.01;

 for(j=0;j<M;j++){
  for(i=0;i<n;i++){
   a[j][i]=s;
  }
 }

 j=0;
 
sp=(spread_CDO(0,3, M,n,r,0.12,0.22,R,a));
  
 do{ 
  j++;
  f=sp-spread[20];
  for(k=26;k<46;k++){
   for(i=0;i<n;i++){
    a[k][i]=s+eps;
    }
  }
  
 sp=(spread_CDO(0,3,M,n,r,0.12,0.22,R,a));
 f1=sp-spread[20];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

 }while(j<Max_iter);
 
 init[4]=s;

nb--;



if(nb>0){

 s=0.1;

 j=0;

for(k=0;k<7;k++){
 for(i=0;i<l1;i++){
  a[k][i]=init[0];
 }
}

for(k=0;k<7;k++){
  for(i=l1;i<n;i++){
   a[k][i]=s;
 }
}


sp=spread_CDO(0,5, M,n,r,0.0,0.03,R,a); 

do{

 f=sp-spread[5];
 j++;

for(k=0;k<7;k++){
 for(i=l1;i<n;i++){
  a[k][i]=s+eps;
 }
}

sp=spread_CDO(0,5,M,n,r,0.0,0.03,R,a);
f1=sp-spread[5];
s1=s-(f/(f1-f))*eps;

 if(fabs(s1-s)<eps) break;

 else s=s1;

}while(j<Max_iter);

init[5]=s1;

s=0.1;

for(k=7;k<13;k++){
   for(i=0;i<l1;i++){
    a[k][i]=init[1];
    }
  }
 

for(j=7;j<13;j++){
 for(i=l1;i<n;i++){
   a[j][i]=s;
  }
 }

j=0;

sp=(spread_CDO(0,5,M,n,r,0.03,0.06,R,a));
  
do{ 
  j++;
  f=sp-spread[9];
  for(k=7;k<13;k++){
   for(i=l1;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,5,M,n,r,0.03,0.06,R,a));
 f1=sp-spread[9];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

 }while(j<Max_iter);

 init[6]=s;

 s=0.01;

for(j=13;j<19;j++){
  for(i=l1;i<n;i++){
   a[j][i]=s;
  }
 }
 
 j=0;
 
 for(k=13;k<19;k++){
   for(i=0;i<l1;i++){
    a[k][i]=init[2];
    }
  }
 
 
sp=(spread_CDO(0,5, M,n,r,0.06,0.09,R,a));
  
do{ 
  j++;
  f=sp-spread[13];
  for(k=13;k<19;k++){
   for(i=l1;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,5,M,n,r,0.06,0.09,R,a));
 f1=sp-spread[13];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

 }while(j<Max_iter);
 
 init[7]=s;
 
 s=0.01;

 for(j=19;j<26;j++){
  for(i=0;i<l1;i++){
    a[j][i]=init[3];
   }
 }

 j=0;

 for(j=19;j<26;j++){
   for(i=l1;i<n;i++){
    a[j][i]=s;
 }
}

 sp=(spread_CDO(0,5, M,n,r,0.09,0.12,R,a));
  
do{ 
  j++;
  f=sp-spread[17];
  for(k=19;k<26;k++){
   for(i=l1;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,5,M,n,r,0.09,0.12,R,a));
 f1=sp-spread[17];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

 }while(j<Max_iter);
 
 init[8]=s;

s=0.01;

for(j=26;j<46;j++){
 for(i=0;i<l1;i++){
   a[j][i]=init[4];
  }
}

for(j=26;j<46;j++){
 for(i=l1;i<n;i++){
   a[j][i]=s;
  }
}

j=0;

 
 sp=(spread_CDO(0,5,M,n,r,0.12,0.22,R,a));
  
 do{ 
  j++;
  f=sp-spread[21];
  for(k=26;k<46;k++){
   for(i=l1;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,5,M,n,r,0.12,0.22,R,a));
 f1=sp-spread[21];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

  }while(j<Max_iter);

  init[9]=s; 
}

nb--;

if(nb>0){

 s=0.1;

 for(k=0;k<7;k++){
  for(i=0;i<l1;i++){
   a[k][i]=init[0];
 }
 for(i=l1;i<l2;i++){
  a[k][i]=init[5];
 }
 for(i=l2;i<n;i++){
  a[k][i]=s;
 }
}

sp=spread_CDO(0,7,M,n,r,0.0,0.03,R,a); 

j=0;

do{

 f=sp-spread[6];
 j++;

for(k=0;k<7;k++){
 for(i=l2;i<n;i++){
  a[k][i]=s+eps;
 }
}
sp=spread_CDO(0,7,M,n,r,0.0,0.03,R,a);
f1=sp-spread[6];
s1=s-(f/(f1-f))*eps;

 if(fabs(s1-s)<eps) break;

 else s=s1;

}while(j<Max_iter);

init[10]=s1;

s=0.1;

for(j=7;j<13;j++){
 for(i=0;i<l1;i++){
   a[j][i]=init[1];
  }
  for(i=l1;i<l2;i++){
   a[j][i]=init[6];
  }
  for(i=l2;i<n;i++){
   a[j][i]=s;
  }
}

j=0;

sp=(spread_CDO(0,7,M,n,r,0.03,0.06,R,a));
  
do{ 
  j++;
  f=sp-spread[10];
  for(k=7;k<13;k++){
   for(i=l2;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,7,M,n,r,0.03,0.06,R,a));
 f1=sp-spread[10];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

 }while(j<Max_iter);

 init[11]=s;

 s=0.01;

for(j=13;j<19;j++){
  for(i=0;i<l1;i++){
   a[j][i]=init[2];
  }
  for(i=l1;i<l2;i++){
   a[j][i]=init[7];
  }
  for(i=l2;i<n;i++){
   a[j][i]=s;
  }
}
 
 j=0;
 
 
sp=(spread_CDO(0,7, M,n,r,0.06,0.09,R,a));
  
do{ 
  j++;
  f=sp-spread[14];
  for(k=13;k<19;k++){
   for(i=l2;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,7,M,n,r,0.06,0.09,R,a));
 f1=sp-spread[14];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

 }while(j<Max_iter);
 
 init[12]=s;
 
 s=0.01;

 for(j=19;j<26;j++){
  for(i=0;i<l1;i++){
    a[j][i]=init[3];
   }
   for(i=l1;i<l2;i++){
    a[j][i]=init[8];
   }
   for(i=l2;i<n;i++){
    a[j][i]=s;
   }
  }

 j=0;

 
 sp=(spread_CDO(0,7, M,n,r,0.09,0.12,R,a));
  
do{ 
  j++;
  f=sp-spread[18];
  for(k=19;k<26;k++){
   for(i=l2;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,7,M,n,r,0.09,0.12,R,a));
 f1=sp-spread[18];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

 }while(j<Max_iter);
 
 init[13]=s;

s=0.01;

for(j=26;j<46;j++){
 for(i=0;i<l1;i++){
   a[j][i]=init[4];
  }
  for(i=l1;i<l2;i++){
   a[j][i]=init[9];
  }
  for(i=l2;i<n;i++){
   a[j][i]=s;
  }
}

j=0;

 sp=(spread_CDO(0,7,M,n,r,0.12,0.22,R,a));
  
 do{ 
  j++;
  f=sp-spread[22];
  for(k=26;k<46;k++){
   for(i=l2;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,7,M,n,r,0.12,0.22,R,a));
 f1=sp-spread[22];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

  }while(j<Max_iter);

  init[14]=s; 
}

nb--;

if(nb>0){

s=0.1;

for(k=0;k<7;k++){
  for(i=0;i<l1;i++){
   a[k][i]=init[0];
 }
 for(i=l1;i<l2;i++){
  a[k][i]=init[5];
 }
 for(i=l2;i<l3;i++){
   a[k][i]=init[10];
  }
  for(i=l3;i<n;i++){
   a[k][i]=s;
  }
}

 sp=spread_CDO(0,10, M,n,r,0.0,0.03,R,a); 

j=0;

do{

 f=sp-spread[7];
 j++;

for(k=0;k<7;k++){
 for(i=l3;i<n;i++){
  a[k][i]=s+eps;
 }
}
sp=spread_CDO(0,10,M,n,r,0.0,0.03,R,a);
f1=sp-spread[7];
s1=s-(f/(f1-f))*eps;

 if(fabs(s1-s)<eps) break;

 else s=s1;

}while(j<Max_iter);

init[15]=s1;

s=0.1;

for(k=7;k<13;k++){
  for(i=0;i<l1;i++){
   a[k][i]=init[1];
 }
 for(i=l1;i<l2;i++){
  a[k][i]=init[6];
 }
 for(i=l2;i<l3;i++){
   a[k][i]=init[11];
  }
  for(i=l3;i<n;i++){
   a[k][i]=s;
  }
}
j=0;

sp=(spread_CDO(0,10,M,n,r,0.03,0.06,R,a));
  
do{ 
  j++;
  f=sp-spread[11];
  for(k=7;k<13;k++){
   for(i=l3;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,10,M,n,r,0.03,0.06,R,a));
 f1=sp-spread[11];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

 }while(j<Max_iter);

 init[16]=s;

 s=0.1;

for(k=13;k<19;k++){
  for(i=0;i<l1;i++){
   a[k][i]=init[2];
 }
 for(i=l1;i<l2;i++){
  a[k][i]=init[7];
 }
 for(i=l2;i<l3;i++){
   a[k][i]=init[12];
  }
  for(i=l3;i<n;i++){
   a[k][i]=s;
  }
}
 
 j=0;
 
 
sp=(spread_CDO(0,10, M,n,r,0.06,0.09,R,a));
  
do{ 
  j++;
  f=sp-spread[15];
  for(k=13;k<19;k++){
   for(i=l3;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,10,M,n,r,0.06,0.09,R,a));
 f1=sp-spread[15];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

 }while(j<Max_iter);
 
 init[17]=s;
 
 s=0.01;

 for(k=19;k<26;k++){
  for(i=0;i<l1;i++){
   a[k][i]=init[3];
 }
 for(i=l1;i<l2;i++){
  a[k][i]=init[8];
 }
 for(i=l2;i<l3;i++){
   a[k][i]=init[13];
 }
  for(i=l3;i<n;i++){
   a[k][i]=s;
  }
}
j=0;

sp=(spread_CDO(0,10, M,n,r,0.09,0.12,R,a));
  
do{ 
  j++;
  f=sp-spread[19];
  for(k=19;k<26;k++){
   for(i=l3;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,10,M,n,r,0.09,0.12,R,a));
 f1=sp-spread[19];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

 }while(j<Max_iter);
 
 init[18]=s;

s=0.01;

for(k=26;k<46;k++){
  for(i=0;i<l1;i++){
   a[k][i]=init[4];
 }
 for(i=l1;i<l2;i++){
  a[k][i]=init[9];
 }
 for(i=l2;i<l3;i++){
   a[k][i]=init[14];
  }
  for(i=l3;i<n;i++){
   a[k][i]=s;
  }
}

j=0;

sp=(spread_CDO(0,10,M,n,r,0.12,0.22,R,a));
  
 do{ 
  j++;
  f=sp-spread[23];
  for(k=26;k<46;k++){
   for(i=l3;i<n;i++){
    a[k][i]=s+eps;
    }
  }
 
 sp=(spread_CDO(0,10,M,n,r,0.12,0.22,R,a));
 f1=sp-spread[23];
 s2=s-(f/(f1-f))*eps;
 
 if(fabs(s2-s)<eps) break;
 
 else s=s2;

  }while(j<Max_iter);

  init[19]=s; 
}


 return init;

 for(i=0;i<M;i++){
  free(a[i]);
 }
 free(a);
 free(init);
}

 





/** this function simulate a path of a transition rate between (0,t)**/

double  ***calib_rates_CDO(double t,double T,int nb,double vol,double r, double* spread){

double** v;
int i,j,d;
double **trans;
double s;
double* sum;
int  L=0;
double ***p;
int nb1;
double a;

double** w;
double ***dym_trans;
double l1,l2,l3;
double* init;
int m;
int M=125;
double R=0.4;

sum=malloc(M*sizeof(double));
w=malloc(M*sizeof(double*));
dym_trans=malloc(M*sizeof(double**));

for(i=0;i<M;i++){
 w[i]=malloc(nb*sizeof(double));
 dym_trans[i]=malloc(nb*sizeof(double*));
}
for(i=0;i<M;i++){
 for(j=0;j<nb;j++){
  dym_trans[i][j]=malloc((nb+1)*sizeof(double));
 }
}


double* drift;
drift=malloc(M*sizeof(double));
v=malloc(M*sizeof(double*));
trans=malloc(M*sizeof(double*));
p=malloc(M*sizeof(double**));




for(j=0;j<M;j++){
 drift[j]=0;
 sum[j]=0;
 v[j]=malloc(M*sizeof(double));
 trans[j]=malloc(nb*sizeof(double));
 p[j]=malloc(M*sizeof(double*));
}

l1=(0.3*nb);
l2=(0.5*nb);
l3=0.7*nb;


if(T<=3){
   nb1=1;
   m=l1;
  }
  else{
     if(T<=5){
        nb1=2;
         m=l2;
      }
     else{
          if(T<=7){
             nb1=3;
             m=l3;
           }
          else{
          nb1=4;
          m=nb;
          }
    }
 }

/** Initialisation des transitions rates **/

init=malloc((5*nb1)*sizeof(double));

init=initial_rates_CDO(T,nb,r,spread);


for(i=0;i<M;i++){
  for(j=0;j<nb;j++){
   trans[i][j]=0.1;
  }
}

if(nb1==1){
 for(j=0;j<nb;j++){
   for(i=0;i<7;i++){
    trans[i][j]=init[0];
   }
   for(i=7;i<13;i++){
    trans[i][j]=init[1];
   }
   for(i=13;i<19;i++){
    trans[i][j]=init[2];
   }

   for(i=19;i<26;i++){
    trans[i][j]=init[3];
   }
   for(i=26;i<46;i++){
    trans[i][j]=init[4];
  }
 }
}
 

if(nb1==2){
 
   for(i=0;i<7;i++){
     for(j=0;j<l1;j++){
       trans[i][j]=init[0];
     }
     for(j=l1;j<nb;j++){
       trans[i][j]=init[5];
     }
   }

   for(i=7;i<13;i++){
     for(j=0;j<l1;j++){
      trans[i][j]=init[1];
    }
    for(j=l1;j<nb;j++){
       trans[i][j]=init[6];
    }
   }

   for(i=13;i<19;i++){
     for(j=0;j<l1;j++){
      trans[i][j]=init[2];
     }

   for(j=l1;j<nb;j++){
       trans[i][j]=init[7];
    }
   }


   for(i=19;i<26;i++){
     for(j=0;j<l1;j++){
       trans[i][j]=init[3];
     }
    for(j=l1;j<nb;j++){
       trans[i][j]=init[8];
    }
   }

   for(i=26;i<46;i++){
     for(j=0;j<l1;j++){
      trans[i][j]=init[4];
    }
    for(j=l1;j<nb;j++){
       trans[i][j]=init[9];
    }
   }
 }

if(nb1==3){

 for(i=0;i<7;i++){
     for(j=0;j<l1;j++){
       trans[i][j]=init[0];
     }
     for(j=l1;j<l2;j++){
       trans[i][j]=init[5];
     }
     for(j=l2;j<nb;j++){
       trans[i][j]=init[10];
     }  
 }

   for(i=7;i<13;i++){
     for(j=0;j<l1;j++){
      trans[i][j]=init[1];
    }
    for(j=l1;j<l2;j++){
       trans[i][j]=init[6];
    }
    for(j=l2;j<nb;j++){
       trans[i][j]=init[11];
    }
 }

   for(i=13;i<19;i++){
     for(j=0;j<l1;j++){
      trans[i][j]=init[2];
     }

   for(j=l1;j<l2;j++){
       trans[i][j]=init[7];
    }
    for(j=l2;j<nb;j++){
       trans[i][j]=init[12];
    }
   }
  for(i=19;i<26;i++){
     for(j=0;j<l1;j++){
       trans[i][j]=init[3];
     }
    for(j=l1;j<l2;j++){
       trans[i][j]=init[8];
    }
    for(j=l2;j<nb;j++){
       trans[i][j]=init[13];
    }
  }

   for(i=26;i<46;i++){
     for(j=0;j<l1;j++){
      trans[i][j]=init[4];
    }
    for(j=l1;j<l2;j++){
       trans[i][j]=init[9];
    }
    for(j=l2;j<nb;j++){
       trans[i][j]=init[14];
    }
  }
}

if(nb1==4){

 for(i=0;i<7;i++){
     for(j=0;j<l1;j++){
       trans[i][j]=init[0];
     }
     for(j=l1;j<l2;j++){
       trans[i][j]=init[5];
     }
     for(j=l2;j<l3;j++){
       trans[i][j]=init[10];
     }  
     for(j=l3;j<nb;j++){
       trans[i][j]=init[15];
     } 
  }

   for(i=7;i<13;i++){
     for(j=0;j<l1;j++){
      trans[i][j]=init[1];
    }
    for(j=l1;j<l2;j++){
       trans[i][j]=init[6];
    }
    for(j=l2;j<l3;j++){
       trans[i][j]=init[11];
    }
    for(j=l3;j<nb;j++){
       trans[i][j]=init[16];
    }
 }

   for(i=13;i<19;i++){
     for(j=0;j<l1;j++){
      trans[i][j]=init[2];
     }

   for(j=l1;j<l2;j++){
       trans[i][j]=init[7];
    }
    for(j=l2;j<l3;j++){
       trans[i][j]=init[12];
    }
    for(j=l3;j<nb;j++){
       trans[i][j]=init[17];
    }
  }
  for(i=19;i<26;i++){
     for(j=0;j<l1;j++){
       trans[i][j]=init[3];
     }
    for(j=l1;j<l2;j++){
       trans[i][j]=init[8];
    }
    for(j=l2;j<l3;j++){
       trans[i][j]=init[13];
    }
    for(j=l3;j<nb;j++){
       trans[i][j]=init[18];
    }
  }
  for(i=26;i<46;i++){
     for(j=0;j<l1;j++){
      trans[i][j]=init[4];
    }
    for(j=l1;j<l2;j++){
       trans[i][j]=init[9];
    }
    for(j=l2;j<l3;j++){
       trans[i][j]=init[14];
    }
    for(j=l3;j<nb;j++){
       trans[i][j]=init[19];
    }
  }
}



for(i=0;i<M;i++){
 for(j=0;j<M;j++){
 p[i][j]=malloc(nb*sizeof(double));
 }
}

/** Initialisation des  probas de passer de i à l'instant t de j à l'instant t+d   **/

int jcompt=0;

for(i=0;i<M;i++){
 for(j=0;j<nb;j++){
   w[i][j]=0;
  }  

 }


int c=L;

p=proba(jcompt*t/nb,T,M,nb,trans);

 do{
if(L!=c) p=proba(jcompt*t/nb,T,M,nb,trans);

/**initialisation des volatilités de transition **/

for(i=0;i<M;i++){
  for(j=0;j<M;j++){
   if(i==j) v[i][j]=-(T-jcompt*t/nb)*vol;
   else     v[i][j]=0;
  }
}


/** initialisation de la perte **/



 a=unif();
 c=L;
 if(a<trans[L][jcompt]*t/nb) L+=1;



/**Initialisation des poids w pour calculer recursivement les vols des probas**/


for(j=L+1;j<M;j++){
 sum[j]=0;
 for(d=0;d<nb;d++){
   s=0;
   for(i=d;i<nb;i++){
    s+=trans[j][i];
   }
   
   s=exp(-((T-d*T/nb)/nb)*s);
   w[j][d]=p[L][j][d]*trans[j][d]*s;
   sum[j]+=w[j][d];  
  }
  if(sum[j]==0) sum[j]=0.00000001;
}

for(j=L+1;j<M;j++){
  for(d=0;d<nb;d++){
    w[j][d]=w[j][d]*1./sum[j];
  }
}

for(j=L+1;j<M;j++){
  v[L][j]=0;
  
  for(d=0;d<nb;d++){
   if(p[L][j-1][nb-1]==0) p[L][j-1][nb-1]=0.000000001;
   v[L][j]+=T/nb*w[j][d]*(v[L][j-1]/p[L][j-1][nb-1]+vol/trans[j-1][d]-vol*(T-d*T/nb));
  }
   v[L][j]=p[L][j][nb-1]*v[L][j];
 }

/** evaluation des drifts des transitions rates **/


for(j=L;j<M;j++){
 if(p[L][j][nb-1]==0) p[L][j][nb-1]=0.00000001;
 drift[j]=-vol*T*v[L][j]*1./p[L][j][nb-1];
}

for(j=0;j<M;j++){
  for(i=0;i<nb;i++){
   dym_trans[j][i][jcompt]=trans[j][i];
   trans[j][i]+=drift[j]*t/nb+vol*(T-T*i/nb)*gaussian()*sqrt(t/nb);
   dym_trans[j][jcompt][nb]=L;
  }
 }
 jcompt++;

}while(jcompt<nb);

 return dym_trans;
}

double  ***calib_rates_CDS(double t,double T,int nb,double vol,double r, double* spread){
int M=125;
double R=0.4;
double** v;
int i,j;
double **trans;

double* sum;

int  L=0;
int k;
double ***p;


int nb1;
double *init_rates;
double** w;
double ***dym_trans;
double l1,l2,l3;

int m;

sum=malloc(M*sizeof(double));
w=malloc(M*sizeof(double*));
dym_trans=malloc(M*sizeof(double**));

for(i=0;i<M;i++){
 w[i]=malloc(nb*sizeof(double));
 dym_trans[i]=malloc(nb*sizeof(double*));
}
for(i=0;i<M;i++){
 for(j=0;j<nb;j++){
  dym_trans[i][j]=malloc((nb+1)*sizeof(double));
 }
}


double* drift;
drift=malloc(M*sizeof(double));
v=malloc(M*sizeof(double*));
trans=malloc(M*sizeof(double*));
p=malloc(M*sizeof(double**));




for(j=0;j<M;j++){
 drift[j]=0;
 sum[j]=0;
 v[j]=malloc(M*sizeof(double));
 trans[j]=malloc(nb*sizeof(double));
 p[j]=malloc(M*sizeof(double*));
}

l1=(0.3*nb);
l2=(0.5*nb);
l3=0.7*nb;


if(T<=3){
   nb1=1;
   m=l1;
  }
  else{
     if(T<=5){
        nb1=2;
         m=l2;
      }
     else{
          if(T<=7){
             nb1=3;
             m=l3;
           }
          else{
          nb1=4;
          m=nb;
          }
    }
 }

/** Initialisation des transitions rates **/

init_rates=malloc((nb1)*sizeof(double));

init_rates=initial_rates_CDS(T,nb,r,spread);

/** Initialisation des transitions rates **/

if(nb1==1){
 for(i=0;i<M;i++){
  for(k=0;k<nb;k++){
  trans[i][k]=init_rates[0];
  }
 }
}

if(nb1==2){
for(i=0;i<M;i++){ 
 for(k=0;k<l1;k++){
  trans[i][k]=init_rates[0];
  }


  for(k=l1;k<nb;k++){
   trans[i][k]=init_rates[1];
   }
  }
 }

if(nb1==3){
 for(i=0;i<M;i++){
  for(k=0;k<l1;k++){
   trans[i][k]=init_rates[0];
  }
  for(k=l1;k<l2;k++){
   trans[i][k]=init_rates[1];
  }
  
 for(k=l2;k<nb;k++){
   trans[i][k]=init_rates[2];
  }
 }
}
 if(nb1==4){
 for(i=0;i<M;i++){
  for(k=0;k<l1;k++){
  trans[i][k]=init_rates[0];
  }
  
  for(k=l1;k<l2;k++){
   trans[i][k]=init_rates[1];
 }
 for(k=l2;k<l3;k++){
   trans[i][k]=init_rates[2];
 }
 for(k=l3;k<nb;k++){
  trans[i][k]=init_rates[3];
  }
 }
}

for(i=0;i<M;i++){
 for(j=0;j<M;j++){
 p[i][j]=malloc(nb*sizeof(double));
 }
}

/** Initialisation des  probas de passer de i à l'instant t de j à l'instant t+d   **/

int jcompt=0;

for(i=0;i<M;i++){
 for(j=0;j<nb;j++){
   w[i][j]=0;
  }  

 }




int c=L;

p=proba(jcompt*t/nb,T,M,nb,trans);

 do{
if(L!=c) p=proba(jcompt*t/nb,T,M,nb,trans);

/**initialisation des volatilités de transition **/

for(i=0;i<M;i++){
  for(j=0;j<M;j++){
   if(i==j) v[i][j]=-(T-jcompt*t/nb)*vol;
   else     v[i][j]=0;
  }
}


/** initialisation de la perte **/

double a;
int d;
double s;
 a=unif();
 c=L;
 if(a<trans[L][jcompt]*t/nb) L+=1;



/**Initialisation des poids w pour calculer recursivement les vols des probas**/


for(j=L+1;j<M;j++){
 sum[j]=0;
 for(d=0;d<nb;d++){
   s=0;
   for(i=d;i<nb;i++){
    s+=trans[j][i];
   }
   
   s=exp(-((T-d*T/nb)/nb)*s);
   w[j][d]=p[L][j][d]*trans[j][d]*s;
   sum[j]+=w[j][d];  
  }
  if(sum[j]==0) sum[j]=0.00000001;
}

for(j=L+1;j<M;j++){
  for(d=0;d<nb;d++){
    w[j][d]=w[j][d]*1./sum[j];
  }
}

for(j=L+1;j<M;j++){
  v[L][j]=0;
  
  for(d=0;d<nb;d++){
   if(p[L][j-1][nb-1]==0) p[L][j-1][nb-1]=0.000000001;
   v[L][j]+=T/nb*w[j][d]*(v[L][j-1]/p[L][j-1][nb-1]+vol/trans[j-1][d]-vol*(T-d*T/nb));
  }
   v[L][j]=p[L][j][nb-1]*v[L][j];
 }

/** evaluation des drifts des transitions rates **/


for(j=L;j<M;j++){
 if(p[L][j][nb-1]==0) p[L][j][nb-1]=0.00000001;
 drift[j]=-vol*T*v[L][j]*1./p[L][j][nb-1];
}

for(j=0;j<M;j++){
  for(i=0;i<nb;i++){
   dym_trans[j][i][jcompt]=trans[j][i];
   trans[j][i]+=drift[j]*t/nb+vol*(T-T*i/nb)*gaussian()*sqrt(t/nb);
   dym_trans[j][jcompt][nb]=L;
  }
 }
 jcompt++;

}while(jcompt<nb);

 return dym_trans;
}
