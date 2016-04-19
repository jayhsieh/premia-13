#include <stdlib.h>
#include<stdio.h>
#include<math.h>
#define M_Pi 3.14159265358979323846264338327950288


#define  N 4
#define nu 20000
#define deltat 0.01
#define horizon 2
#define TMAX  10000000
double sup[TMAX];
double P[TMAX];


 double fplus(double S, double K)
{
  if(K>S){
	return (K-S);
  } else {
	return 0;
  }
}

 double drand48(void);



 double repartition(double x)
{
  double p, b1, b2, b3, b4, b5, t, w;
 
  p=0.2316419;
  b1=0.319381530;
  b2=-0.356563782;
  b3=1.781477937;
  b4=-1.821255978;
  b5=1.330274429;
  t=1/(1+p*fabs(x));
  w=1-exp(-0.5*x*x)*(b1*t+b2*pow(t,2)+b3*pow(t,3)+b4*pow(t,4)+b5*pow(t,5))/sqrt(2*M_Pi);
  
 if (x>0) return w;
 else return 1-w;
 
}

 double putBlackScholes(double x, double spot,double K,double r,double horizon1)
{
  double mu=(r+x*x/2);
  double sigmasqrt=x*sqrt(horizon1);
  double d1=(log(spot/K)+mu*horizon1)/(sigmasqrt);
  double d2=d1-sigmasqrt;


  
 return -spot*repartition(-d1)+K*exp(-r*horizon1)*repartition(-d2);
 
}

	
double*  tiragedetrajectoires(double spot,double r)
{
  
  
   int i;
   int t=0;
   int n=(int)floor(horizon/deltat);
   double pasdetemps=sqrt(deltat);
   double alpha=0.1;
   double beta=1;
   double rho=0.5;
   double U1;
   double U2;
   double g1;
   double g2;
   double volatilite;
   double sigma0=0.2;
   double accroissement1;
   double accroissement2;
   double cours;
   double *tirages=malloc((1+n*nu)*sizeof(double));
   
   

   for (i=0;i<nu;i++)
{
  /*simulation de deux gaussiennes et des accroissements browniens correles correspondant*/

 
  cours=spot;
  volatilite=0.2;
	/*	printf("accroissement1\n");
		printf("%g\n", accroissement1);*/
  
	
/*simulation du cours*/
/*simulation de la volatilité*/	
 for (t=0;t<n;t++){
     U1=drand48();
	 U2=drand48();
	 g1=sqrt(-2*log(U1))*cos(2*M_Pi*U2);
	 g2=sqrt(-2*log(U1))*sin(2*M_Pi*U2);
	 accroissement1=pasdetemps*g1;
	 accroissement2=pasdetemps*(rho*g1+sqrt(1-rho*rho)*g2);
	 cours=cours*((1+deltat*r)+volatilite*accroissement1);
	 volatilite=sigma0+(volatilite-sigma0)*(1-deltat*beta)+alpha*accroissement2;
	 tirages[t+i*n]=cours;
	 /* printf("%g\n",tirages[t+i*n]);*/
	 
	 
 }
}
   return tirages;
	
}

double* donnees(double*tirages,double spot, double r)
{

   int n=(int)floor(horizon/deltat);
   int i;
   double sigma0=0.3;
   int max=floor(((1.3*spot)-(0.7*spot))/5);
   int j;
   int k;
   int m;
   double horizon1;
   int l;
   double K=0.7*spot;
   double *G=malloc((1+nu*N*max)*sizeof(double));
   
   
   
   for (m=0;m<max;m++){
	 
	 K=K+5;
	 horizon1=0;
	

	 for (j=0;j<N;j++){
	   
	   horizon1=horizon1+0.25;
	   l=(int)floor(horizon1/deltat);
	   P[j+m*N]=putBlackScholes(sigma0,spot, K, r, horizon1);
	   /* printf("ok\n");
		  printf("%g %g %g\n",K,horizon1,P[j+m*N]);*/
	   for (i=0;i<nu;i++){
		 G[j+m*N+i*max*N]=exp(-r*horizon1)*fplus(tirages[l+i*n],K)-P[j+m*N];
		 /* printf("%g\n",G[j+m*N+i*max*N]);*/
		 

		 if(fabs(G[j+m*N+i*max*N])>sup[j+m*N]){

		   sup[j+m*N]= fabs(G[j+m*N+i*max*N]);
		 }
		 /* printf("%g\n",sup[j+m*N]);*/
		 
	   }

	   for(k=0;k<nu;k++){
		 if(sup[j+m*N]>0){
		   G[j+m*N+k*max*N]=(1/sup[j+m*N])* G[j+m*N+k*max*N];
		 } 
		 /* printf("%g %g\n",sup[j+m*N],G[j+m*N+k*max*N]);*/
	   }
	 }
	 }
   
   return G;
   
}






double fonctionnelle(double *G,double *tirages,double *lambda, double spot, double r)
{
  double S;
  double W=0;
  int i;
  int j;
  int m;
  double inverse=(double)1/nu;
  int max=floor(((1.3*spot)-(0.7*spot))/5);
  
  

  for(i=0;i<nu;i++){S=0;

  for (m=0;m<max;m++){

	for(j=0;j<N;j++){
	  /* printf("%d %d %g %g\n",j+m*N,j+m*N+max*N*i,lambda[j+m*N],G[j+m*N+max*N*i]);*/
		S=S+lambda[j+m*N]*G[j+m*N+max*N*i];
	
	  
	}
  }
  
	  
	W=W+inverse*exp(S);
	/*	printf("%g\n",W);*/
	
  }
  
  return log(W);
}


int gradient(double *G,double *Grad,double *tirages,double*lambda,double spot, double r)
{
  int i;
  int j;
  int k;
  int l;
  int m;
  double S;
  double Z;
  double resultat;
  int max=floor(((1.3*spot)-(0.7*spot))/5);
  

  for(k=0;k<max;k++){
	for(l=0;l<N;l++){
	Z=0;
	resultat=0;
	for(i=0;i<nu;i++){
	S=0;
	for (m=0;m<max;m++){

	  for(j=0;j<N;j++){	
		
		S=S+lambda[j+m*N]*G[j+m*N+max*N*i];
	  }
	}
	Z=Z+exp(S);
	resultat=resultat+G[l+k*N+max*N*i]*exp(S);
	}
	Grad[l+k*N]=(1/Z)*resultat;
	}
  }
  
  return 0;
  }

double normeinfinie(double *x,double *y,double spot,double r)
{
  int i;
  double alpha=0;
  int max=floor(((1.3*spot)-(0.7*spot))/5);

  for(i=0;i<N*max;i++){
	if(fabs(x[i]-y[i])>alpha){
	  alpha=fabs(x[i]-y[i]);
	}
  }
  
  return alpha;
}


double produitscalaire(double *x,double *y,double spot,double r)
{
  int i;
  double alpha=0;
  int max=floor(((1.3*spot)-(0.7*spot))/5);
  

  for(i=0;i<N*max;i++){
	
	alpha=alpha+x[i]*y[i];
	
  }
  return alpha;
}

int addition(double*add,double *x,double *y,double spot,double r)
{ int i;

  int max=floor(((1.3*spot)-(0.7*spot))/5);

  for(i=0;i<N*max;i++){
	
	add[i]=x[i]+y[i];
  }
  
  return 0;
  
}

int multiplication(double*mult,double a,double *x,double spot, double r)
{
  int i;
  int max=floor(((1.3*spot)-(0.7*spot))/5);
  

  for(i=0;i<N*max;i++){
	mult[i]=a*x[i];
	
  }
  return 0;
  
}

double fonctiondemerite(double *x,double t,double *d,double *G,double *tirages,double spot,double r)
{ int max=floor(((1.3*spot)-(0.7*spot))/5);
  double *add=malloc((N*max+1)*sizeof(double));
  double *mult=malloc((N*max+1)*sizeof(double));
  multiplication(mult,-t,d,spot,r);
  addition(add,x,mult,spot,r);


  return fonctionnelle(G,tirages,add, spot, r);
}

  

double  RegledeWolfe(double *x,double*d,double t,double*G,double *tirages,double spot,double r)
	 
{ int max=floor(((1.3*spot)-(0.7*spot))/5);
  double m1=0.2;
  double m2=0.7;
  double td=10;
  double tg=0;
  double a=2;
  double *Grad=malloc((N*max+1)*sizeof(double));
  double *grad=malloc((N*max+1)*sizeof(double));
  double *add=malloc((N*max+1)*sizeof(double));
  double *mult=malloc((N*max+1)*sizeof(double));
  multiplication(mult,-t,d,spot,r);
  addition(add,x,mult,spot,r);
  gradient(G,Grad,tirages,x ,spot,r);
  gradient(G,grad,tirages,add ,spot,r);
	
  while(((fonctiondemerite(x,t,d,G,tirages,spot,r)>fonctiondemerite(x,0,d,G,tirages,spot,r)-m1*t*produitscalaire(Grad,Grad,spot,r))&&(fabs(td-tg)>1E-6))||((produitscalaire(Grad,grad,spot,r)>m2*produitscalaire(Grad,Grad,spot,r))&&(fabs(td-tg)>1E-6))){printf("no\n");printf("%f %f\n",tg,td);
  
  	while (td==0){
	  printf("ok\n");
	  if (fonctiondemerite(x,t,d,G,tirages,spot,r)>fonctiondemerite(x,0,d,G,tirages,spot,r)-m1*t*produitscalaire(Grad,Grad,spot,r)){td=t;
	} else {
	t=a*t;
	}
	}

	t=0.5*(td+tg);
	
	
	if(fonctiondemerite(x,t,d,G,tirages,spot,r)>fonctiondemerite(x,0,d,G,tirages,spot,r)-m1*t* produitscalaire(Grad,Grad,spot,r)){
	td=t;
	t=0.5*(tg+td);
	} else {
	  if((fonctiondemerite(x,t,d,G,tirages,spot,r)<=fonctiondemerite(x,0,d,G,tirages,spot,r)-m1*t* produitscalaire(Grad,Grad,spot,r))&&( produitscalaire(Grad,grad,spot,r)>m2* produitscalaire(Grad,Grad,spot,r))){
	
  	tg=t;
	t=0.5*(tg+td);
	  }
	}

	multiplication(mult,-t,d,spot,r);
	addition(add,x,mult,spot,r);
	gradient(G,grad,tirages,add ,spot,r);

	}
  
		 
  free(Grad);
  free(grad);
  free(add);
  free(mult);
  
  
return t;}



int methodedugradient(double *G,double *lambda,double *tirages,double spot, double r )
{
  int i;
  double pas;
  double test;
  double t;
  int max=floor(((1.3*spot)-(0.7*spot))/5);
  double *Grad=malloc((N*max+1)*sizeof(double));
  double *grad=malloc((N*max+1)*sizeof(double));
  double *mult=malloc((N*max+1)*sizeof(double));
  double *add=malloc((N*max+1)*sizeof(double));
  double *Mult=malloc((N*max+1)*sizeof(double));
  double *lambda0=malloc((N*max+1)*sizeof(double));
  double*zero=malloc((N*max+1)*sizeof(double));
  double *direction=malloc((N*max+1)*sizeof(double));
  for(i=0;i<N*max;i++)
	{
	  lambda0[i]=1;
	  zero[i]=0;
	} 
	  multiplication(lambda,0,lambda0,spot,r);
	  gradient(G,Grad,tirages,lambda ,spot,r);
	  gradient(G,direction,tirages,lambda ,spot,r);
	  test=normeinfinie(Grad,lambda,spot,r);


	  while ((normeinfinie(lambda0,lambda,spot,r)>1E-10)&&/*(fabs(pasapriori*produitscalaire(Grad,Grad,spot,r)-fonctionnelle(G,tirages,lambda,spot,r)+fonctionnelle(G,tirages,add,spot,r))>1E-20)*/(normeinfinie(zero,Grad,spot,r)/test>1E-5)){
		printf("W\n");
		t=1;
		pas=RegledeWolfe(lambda,direction,t,G,tirages,spot,r);
		
	
	for(i=0;i<N*max;i++){
	  lambda0[i]=lambda[i];
	
  }
	multiplication(Mult,- pas,direction,spot,r);
	addition(lambda,lambda0,Mult,spot,r);
	gradient(G,Grad,tirages,lambda ,spot,r);
	gradient(G,grad,tirages,lambda0 ,spot,r);
	multiplication(mult,produitscalaire(Grad,Grad,spot,r)/produitscalaire(grad,grad,spot,r),Grad,spot,r);
	addition(direction,direction,mult,spot,r);
	printf("%g\n",fonctionnelle(G,tirages, lambda, spot,r));
	  }
	  

	  free(Grad);
	  free(lambda0);
	  free(add);
	  free(mult);
	  free(Mult);
	  free(zero);
	  free(direction);
	  free(grad);
	
  return 0;
}


int poids(double *G,double*poids,double *tirages,double spot,double r)
{
  int i;
  int j;
  int m;
  double S;
  double Z=0;
  double somme=0;
  int max=floor(((1.3*spot)-(0.7*spot))/5);
  double *l=malloc((1+N*max)*sizeof(double));
  methodedugradient(G,l,tirages,spot,r);
  
  /*  printf("Z\n");
	  printf("%g\n",Z);*/
  
  for(i=0;i<nu;i++){
	S=0;
	for (m=0;m<max;m++){

	  for(j=0;j<N;j++){	/*printf("%g %g\n",l[j+m*N],G[j+m*N+max*N*i]);*/
		
		S=S+l[j+m*N]*G[j+m*N+max*N*i];
	  }
	}
	poids[i]=/*(1/Zdelambda(lambda0,lambda,spot,r))**/exp(S);
	Z=Z+poids[i];
	
  }
  for (i=0;i<nu;i++){
	poids[i]=(1/Z)*poids[i];
	somme=somme+poids[i];
	printf("%d %g\n",i,poids[i]);
  }
  
  free(l);
  
  printf("somme\n");
  printf("%g\n",somme);
  return 0;
  
}

  



int main(int argc,char**argv)
{
  double spot=atof(argv[1]);
  double r=atof(argv[2]);
  int i;
  int j;
  int m;
  int n=(int)floor(horizon/deltat);
  double W;
  int max=floor(((1.3*spot)-(0.7*spot))/5);
  double *tirages=malloc((1+nu*n)*sizeof(double));
  /* double *lambda1=malloc((N*max+1)*sizeof(double));*/
  double *proba=malloc((1+nu)*sizeof(double));
  double *G=malloc((1+nu*max*N)*sizeof(double));
  /*double *resultat=malloc((1+nu*max*N)*sizeof(double));*/

  tirages=tiragedetrajectoires(spot,r);
	  /* W=fonctionnelle(lambda, spot, r);
		 printf("%g\n",W);*/

  /* lambda1=methodedugradient(tirages,spot,r);*/
   G=donnees(tirages,spot,r);
   poids(G,proba,tirages,spot,r);

   
		
		for (m=0;m<max;m++){
		  for(j=0;j<N;j++){
			 W=0;
			 for(i=0;i<nu;i++){
			  
			W=W+proba[i]*G[j+m*N+max*N*i];
	
	  
	}
		if(fabs(W)>1E-3){
		  printf("%d %g\n",j+m*N,W);
		}
		  }
		  
		
		}

  free(proba);
  free(G);
   /* free(lambda1);*/

  free(tirages);
	  return 0;
}
