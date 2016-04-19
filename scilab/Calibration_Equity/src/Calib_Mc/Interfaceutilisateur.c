#include <stdlib.h>
#include<stdio.h>
#include<math.h>
#define M_Pi 3.14159265358979323846264338327950288

#define nu 20000
#define TMAX  10000000
/*double P[TMAX];
double K[TMAX];
double T[TMAX];*/
double sup[TMAX];


/*This function computes (K-S_T)+  */

 double fplus(double S, double K)
{
  if(K>S){
	return (K-S);
  } else {
	return 0;
  }
}

 double drand48(void);


/* Repartition function of N2(0,Id) */

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


/* Computes the Black Schole's put formula */
 double putBlackScholes(double x, double spot,double K,double r,double horizon1)
{
  double mu=(r+x*x/2);
  double sigmasqrt=x*sqrt(horizon1);
  double d1=(log(spot/K)+mu*horizon1)/(sigmasqrt);
  double d2=d1-sigmasqrt;


  
 return -spot*repartition(-d1)+K*exp(-r*horizon1)*repartition(-d2);
 
}

/* Simulation of the paths */	
double*  tiragedetrajectoires(double alpha,double beta, double rho,double sigma0,double spot,double r,double q,int N,double deltat,double horizon)
{
  
  
   int i;
   int t=0;
   int n=(int)floor(horizon/deltat);
   double pasdetemps=sqrt(deltat);
   double U1;
   double U2;
   double g1;
   double g2;
   double volatilite;
   double accroissement1;
   double accroissement2;
   double cours;
   double *tirages=malloc((1+n*nu)*sizeof(double));
   
   

   for (i=0;i<nu;i++)
{
  /*simulation de deux gaussiennes et des accroissements browniens correles correspondant*/

 
  cours=spot;
  volatilite=sigma0;
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
	 cours=cours*((1+deltat*(r-q))+volatilite*accroissement1);
	 volatilite=sigma0+(volatilite-sigma0)*(1-deltat*beta)+alpha*accroissement2;
	 tirages[t+i*n]=cours;
	 /* printf("%g\n",tirages[t+i*n]);*/
	 
	 
 }
}
   return tirages;
	
}


/*Renormalization of the payoffs */

double* donnees(double*tirages,double spot,double sigma0, double r,int N,double deltat,double horizon,double *T, double *P, double*K)
{

   int n=(int)floor(horizon/deltat);
   int i;
   int k;
   int m;
   double horizon1;
   int l;
   double *G=malloc((1+nu*N)*sizeof(double));
   
   
   
   for (m=0;m<N;m++){
	 
	 horizon1=T[m];
	 l=(int)floor(horizon1/deltat);
	 
	   for (i=0;i<nu;i++){
		 G[m+i*N]=exp(-r*horizon1)*fplus(tirages[l+i*n],K[m])-P[m];
		

		 if(fabs(G[m+i*N])>sup[m]){

		   sup[m]= fabs(G[m+i*N]);
		 }
		
	   }

	   for(k=0;k<nu;k++){
		 if(sup[m]>0){
		   G[m+k*N]=(1/sup[m])* G[m+k*N];
		 } 
		
	   }
	 }
   
   return G;
   
}





/* computes sup(|x[i]-y[i]);*/

double normeinfinie(double *x,double *y,double spot,double r,int N)
{
  int i;
  double alpha=0;
 

  for(i=0;i<N;i++){
	if(fabs(x[i]-y[i])>alpha){
	  alpha=fabs(x[i]-y[i]);
	}
  }
  
  return alpha;
}

/* computes x.y */
double produitscalaire(double *x,double *y,double spot,double r,int N)
{
  int i;
  double alpha=0;
 

  for(i=0;i<N;i++){
	
	alpha=alpha+x[i]*y[i];
	
  }
  return alpha;
}

/* computes x+y */
int addition(double*add,double *x,double *y,double spot,double r,int N)
{ int i;

  

  for(i=0;i<N;i++){
	
	add[i]=x[i]+y[i];
  }
  
  return 0;
  
}

/* computes a*x where a is a real */

int multiplication(double*mult,double a,double *x,double spot, double r,int N)
{
  int i;
 

  for(i=0;i<N;i++){
	mult[i]=a*x[i];
	
  }
  return 0;
  
}

/* computes  W(lambda) */
double fonctionnelle(double *w,double *G,double *tirages,double *lambda, double spot, double r,int N)
{
  double S;
  double W=0;
  int i;
  int j;
  double inverse=(double)1/nu;
 
  
  

  for(i=0;i<nu;i++){S=0;

	for(j=0;j<N;j++){
	 
		S=S+lambda[j]*G[j+i*N];
	
	  
	}
	W=W+inverse*exp(S);
  }
  
  W=log(W);
  
	for(j=0;j<N;j++){
	  W=W+0.5*w[j]*lambda[j]*lambda[j]/(sup[j]*sup[j]);
	}
  return W;
}


/* computes grad(W) */
int gradient(double*w,double *G,double *Grad,double *tirages,double*lambda,double spot, double r,int N)
{
  int i;
  int j;
  int k;
  double S;
  double Z;
  double resultat;
 

  for(k=0;k<N;k++){

	Z=0;
	resultat=0;
	for(i=0;i<nu;i++){
	S=0;

	  for(j=0;j<N;j++){	
		
		S=S+lambda[j]*G[j+i*N];
	  }
	


	Z=Z+exp(S);
	resultat=resultat+G[k+N*i]*exp(S);
	}
	Grad[k]=(1/Z)*resultat+w[k]*lambda[k]/(sup[k]*sup[k]);
	}
  
  return 0;
  }


/* computes W(x+t*d) where x,d are vectors and t is a real */

double fonctiondemerite(double*w,double *x,double t,double *d,double *G,double *tirages,double spot,double r,int N)
{ double *add=malloc((N+1)*sizeof(double));
  double *mult=malloc((N+1)*sizeof(double));
  multiplication(mult,-t,d,spot,r,N);
  addition(add,x,mult,spot,r,N);


  return fonctionnelle(w,G,tirages,add, spot, r,N);
}

  
/* computes the "Regle de Wolfe" */
double  RegledeWolfe(double*w,double *x,double*d,double t,double*G,double *tirages,double spot,double r,int N)
	 
{ double m1=0.4;
  double m2=0.8;
  double td=0;
  double tg=0;
  double a=2;
  double *Grad=malloc((N+1)*sizeof(double));
  double *grad=malloc((N+1)*sizeof(double));
  double *add=malloc((N+1)*sizeof(double));
  double *mult=malloc((N+1)*sizeof(double));
  multiplication(mult,-t,d,spot,r,N);
  addition(add,x,mult,spot,r,N);
  gradient(w,G,Grad,tirages,x ,spot,r,N);
  gradient(w,G,grad,tirages,add ,spot,r,N);
	
  while((fonctiondemerite(w,x,t,d,G,tirages,spot,r,N)>fonctiondemerite(w,x,0,d,G,tirages,spot,r,N)-m1*t*produitscalaire(Grad,Grad,spot,r,N))||(produitscalaire(Grad,grad,spot,r,N)>m2*produitscalaire(Grad,Grad,spot,r,N))){printf("no\n");printf("%f %f\n",tg,td);
  
  	while (td==0){
	  printf("ok\n");
	  if (fonctiondemerite(w,x,t,d,G,tirages,spot,r,N)>fonctiondemerite(w,x,0,d,G,tirages,spot,r,N)-m1*t*produitscalaire(Grad,Grad,spot,r,N)){td=t;
	} else {
	t=a*t;
	}
	}

	t=0.5*(td+tg);
	
	
	if(fonctiondemerite(w,x,t,d,G,tirages,spot,r,N)>fonctiondemerite(w,x,0,d,G,tirages,spot,r,N)-m1*t* produitscalaire(Grad,Grad,spot,r,N)){
	td=t;

	} else {
	 
	
  	tg=t;
	
	  }
	if(fabs(td-tg)<1E-6){
	  return tg;
	} else {
	  t=0.5*(td+tg);
	}

	multiplication(mult,-t,d,spot,r,N);
	addition(add,x,mult,spot,r,N);
	gradient(w,G,grad,tirages,add ,spot,r,N);

	}		 
  free(Grad);
  free(grad);
  free(add);
  free(mult);
  
  return t;}


/* computes a conjugate gradient optimisation */
int methodedugradient(double*w,double *G,double *lambda,double *tirages,double spot, double r,int N )
{
  int i;
  double pas;
  double test;
  double t;
  double *Grad=malloc((N+1)*sizeof(double));
  double *grad=malloc((N+1)*sizeof(double));
  double *mult=malloc((N+1)*sizeof(double));
  double *add=malloc((N+1)*sizeof(double));
  double *Mult=malloc((N+1)*sizeof(double));
  double *lambda0=malloc((N+1)*sizeof(double));
  double*zero=malloc((N+1)*sizeof(double));
  double *direction=malloc((N+1)*sizeof(double));
  for(i=0;i<N;i++)
	{
	  lambda0[i]=1;
	  zero[i]=0;
	} 
	  multiplication(lambda,0,lambda0,spot,r,N);
	  gradient(w,G,Grad,tirages,lambda ,spot,r,N);
	  gradient(w,G,direction,tirages,lambda ,spot,r,N);
	  test=normeinfinie(Grad,lambda,spot,r,N);


	  while ((normeinfinie(lambda0,lambda,spot,r,N)>1E-10)&&(normeinfinie(zero,Grad,spot,r,N)/test>1E-5)){
		printf("W\n");
		t=1;
		pas=RegledeWolfe(w,lambda,direction,t,G,tirages,spot,r,N);
		
	
	for(i=0;i<N;i++){
	  lambda0[i]=lambda[i];
	
  }
	multiplication(Mult,- pas,direction,spot,r,N);
	addition(lambda,lambda0,Mult,spot,r,N);
	gradient(w,G,Grad,tirages,lambda ,spot,r,N);
	gradient(w,G,grad,tirages,lambda0 ,spot,r,N);
	multiplication(mult,produitscalaire(Grad,Grad,spot,r,N)/produitscalaire(grad,grad,spot,r,N),Grad,spot,r,N);
	addition(direction,direction,mult,spot,r,N);
	printf("%g\n",fonctionnelle(w,G,tirages, lambda, spot,r,N));
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


int methodedugradient2(double *w,double *G,double *lambda,double *tirages,double spot, double r,int N )
{
  int i;
  double pas;
  double test;
  double pasapriori=0.0001;
  int max=floor(((1.3*spot)-(0.7*spot))/5);
  double *Grad=malloc((N+1)*sizeof(double));
  double *mult=malloc((N+1)*sizeof(double));
  double *add=malloc((N+1)*sizeof(double));
  double *Mult=malloc((N+1)*sizeof(double));
  double *lambda0=malloc((N+1)*sizeof(double));
  double*zero=malloc((N+1)*sizeof(double));

  for(i=0;i<N;i++)
	{
	  lambda0[i]=1;
	  zero[i]=0;
	} 
  /* multiplication(lambda,0,lambda0,spot,r,N);*/
  methodedugradient(w,G,lambda,tirages,spot,r,N);
  gradient(w,G,Grad,tirages,lambda ,spot,r,N);
  multiplication(mult,- pasapriori,Grad,spot,r,N);
  addition(add,lambda,mult,spot,r,N);
  test=normeinfinie(Grad,lambda,spot,r,N);


	  while ((normeinfinie(lambda0,lambda,spot,r,N)>1E-10)&&(normeinfinie(zero,Grad,spot,r,N)/test>1E-5)){
		printf("W2\n");
		if((fonctionnelle(w,G,tirages,lambda,spot,r,N)>fonctionnelle(w,G,tirages,lambda0,spot,r,N))){
		  lambda=lambda0;
		  pasapriori=0.001*pasapriori;}
		else {if(fabs(fonctionnelle(w,G,tirages,lambda,spot,r,N)-fonctionnelle(w,G,tirages,lambda0,spot,r,N))<1E-5){
		  pasapriori=0.001*pasapriori;
		} else {
		  pasapriori=1;
		}
		}		

	pas=pasapriori*produitscalaire(Grad,Grad,spot,r,N)/(2*(pasapriori*produitscalaire(Grad,Grad,spot,r,N)-fonctionnelle(w,G,tirages,lambda,spot,r,N)+fonctionnelle(w,G,tirages,add,spot,r,N)));

	
	for(i=0;i<N*max;i++){
	  lambda0[i]=lambda[i];
	
  }
	
	if(pas>=0){
	  if((fonctionnelle(w,G,tirages,add,spot,r,N)-fonctionnelle(w,G,tirages,lambda,spot,r,N))<-0.5*pasapriori*produitscalaire(Grad,Grad,spot,r,N)){
		pas=1;
		multiplication(Mult,- pas*pasapriori,Grad,spot,r,N);
		addition(lambda,lambda0,Mult,spot,r,N);
	  } else {
		multiplication(Mult,-pas*pasapriori,Grad,spot,r,N);
		printf("ok\n");
		addition(lambda,lambda0,Mult,spot,r,N);
	  }
	  
	  
	} else {
	  
	  printf("non!\n");
	  pas=1;
	  multiplication(Mult,- pas*pasapriori,Grad,spot,r,N);
	  addition(lambda,lambda0,Mult,spot,r,N);
	  
	}
	 gradient(w,G,Grad,tirages,lambda ,spot,r,N);
	 multiplication(mult,- pasapriori,Grad,spot,r,N);
	 addition(add,lambda,mult,spot,r,N);
	printf("%g\n",fonctionnelle(w,G,tirages, lambda, spot,r,N));
	   
   }
	  free(Grad);
	  free(lambda0);
	  free(add);
	  free(mult);
	  free(Mult);
	  free(zero);
	
  return 0;
}


/* computes the weights */
int poids(double*w,double *G,double*poids,double *tirages,double spot,double r,int N)
{
  int i;
  int j;
  double S;
  double Z=0;
  double somme=0;
  double *l=malloc((1+N)*sizeof(double));
  methodedugradient(w,G,l,tirages,spot,r,N);
  
  for(i=0;i<nu;i++){
	S=0;

	  for(j=0;j<N;j++){	/*printf("%g %g\n",l[j+m*N],G[j+m*N+max*N*i]);*/
		
		S=S+l[j]*G[j+N*i];
	  }
	
	poids[i]=exp(S);
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

/* loads the parameters */  
void loadParameters(char *name_in_parameters,double*  spot, double*  r, double*  q, int*  option_volatility, int* option_step,double*  horizon,double*  sigma0,double * deltat, int *N,double*alpha, double *beta,double*rho,int*option_constraint){
 

  FILE *fic_in_parameters;
  fic_in_parameters = fopen(name_in_parameters,"r");
  fscanf(fic_in_parameters," %lf\n",spot);
  printf(" %f\n",spot[0]);
  fscanf(fic_in_parameters," %lf\n",r);
  printf(" %f\n",r[0]);
  fscanf(fic_in_parameters," %lf\n",q);
  printf(" %f\n",q[0]);
  fscanf(fic_in_parameters," %d\n",option_volatility);
  printf(" %d\n",option_volatility[0]);
  fscanf(fic_in_parameters," %d\n",option_step);
  printf(" %d\n",option_step[0]);
  fscanf(fic_in_parameters," %lf\n",horizon);
  printf(" %f\n",horizon[0]);
  fscanf(fic_in_parameters," %lf\n",sigma0);
  printf(" %f\n",sigma0[0]);
  fscanf(fic_in_parameters," %lf\n",deltat);
  printf(" %f\n",deltat[0]);
  fscanf(fic_in_parameters,"%d\n",N);
  printf(" %d\n",N[0]);
  fscanf(fic_in_parameters,"%lf\n",alpha);
  printf("%f\n",alpha[0]);
  fscanf(fic_in_parameters,"%lf\n",beta);
  printf("%f\n",beta[0]);
  fscanf(fic_in_parameters,"%lf\n",rho);
  printf("%f\n",rho[0]);
  fscanf(fic_in_parameters,"%d\n",option_constraint);
  printf("%d\n",option_constraint[0]);
  
		 




 
 
  fclose(fic_in_parameters);
}

/* loads the data prices */
void loadDataPrices(char *name_in_data,double *K,double *P,double *T,int N){

  FILE *fic_in_data;
  int i;
 
  
  fic_in_data = fopen(name_in_data,"r");

/*    for(i=0;i<N;i++){ */

/*  	fscanf(fic_in_data,"%lf %lf %lf\n",&K[i],&T[i],&P[i]); */
/*  	printf("%f %f %f\n",K[i],T[i],P[i]); */
/*    } */

/*     fclose(fic_in_data); */

/*  } */


/*  /* computes the implied volatility */ */
/*  double volatiliteimplicite (double p,double spot,double K,double r,double horizon) */
/*  {  */
/*   double borneinf=0; */
/*   double bornesup=1; */
/*   double resultat=0.5; */
 
/*   while(fabs(p-putBlackScholes(resultat,spot, K, r, horizon))>0.0001){ */
 
/*     if (putBlackScholes(resultat,spot, K,r,horizon)>p) */
/*  	 { */
/*  	  	bornesup=resultat; */
/*  		borneinf=borneinf; */
/*  		resultat=0.5*(bornesup+borneinf); */
	
/*        } */
/*     else */
/*  	  { */
	
/*  		bornesup=bornesup; */
/*  		borneinf=resultat; */
/*  		resultat=0.5*(bornesup+borneinf);  */
	
/*  	  } */
/*   } */
   
/*     return  resultat; */

 
/*  } */


/*  int main(int argc,char**argv) */
/*  { */
/*    double spot; */
/*    double r; */
/*    double q; */
/*    int i; */
/*    int j; */
/*    int n; */
/*    double tmin; */
/*    double sigma0; */
/*    double variance; */
/*    int N; */
/*    double horizon; */
/*    double deltat; */
/*    double W; */
/*    double alpha; */
/*    double rho; */
/*    double beta; */
/*    char * name_in_data; */
/*    char * name_in_parameters; */
/*    int option_volatility; */
/*    int option_step; */
/*    int option_constraint; */
/*    double*G; */
/*    double*K; */
/*    double*T; */
/*    double*P; */
/*    double*w; */
/*    double *tirages=malloc((1+nu*n)*sizeof(double)); */
/*    double *proba=malloc((1+nu)*sizeof(double)); */
 

/*    name_in_data = "datas.in"; */
/*    name_in_parameters = "parameters.in"; */
/*    loadParameters(name_in_parameters,&spot, &r, &q, &option_volatility, &option_step,&horizon,&sigma0,&deltat,&N,&alpha,&beta,&rho,&option_constraint ); */
 
/*    N=20; */
/*    G=malloc((1+nu*N)*sizeof(double)); */
/*    K=malloc((1+N)*sizeof(double)); */
/*    T=malloc((1+N)*sizeof(double)); */
/*    P=malloc((1+N)*sizeof(double)); */
/*    w=malloc((N+1)*sizeof(double)); */
/*    loadDataPrices(name_in_data,K,P,T,N); */

/*    for(i=0;i<N;i++){if(i==0){ */
/*  	tmin=T[0]; */
/*    } else { */
	
/*  	if(tmin>T[i]){ */
/*  	  tmin=T[i]; */
/*  	} */
/*    } */
/*    } */
/*    printf("tmin\n"); */
/*    printf("%f\n",tmin); */
  
  


  
/*     n=(int)floor(horizon/deltat); */
   
/*     if(option_volatility==1){ */
/*  	 sigma0=0; */
/*  	 variance=0; */
/*  	 rho=0.7; */
/*  	 for(i=0;i<N;i++){ */
/*  	   printf("%d vol\n",i */);
	   printf("%f\n",volatiliteimplicite(P[i],spot,K[i],r,T[i]));
	   sigma0=sigma0+volatiliteimplicite(P[i],spot,K[i],r,T[i]);
	 variance=variance+volatiliteimplicite(P[i],spot,K[i],r,T[i])*volatiliteimplicite(P[i],spot,K[i],r,T[i]);
	
	 }
	 sigma0=(1/(double) N)*sigma0;
	 variance=(1/(double) N)*variance-sigma0*sigma0;
	 alpha=1;
	 beta=2*sqrt(variance)/(1-exp(-tmin));
	 printf("beta\n");
	 printf("%f\n",beta);
	 printf("sigma0\n");
	 printf("%f\n",sigma0);
	
   }
   if(option_step==1){
	 deltat=0.01;
   }
   
   if(option_constraint==1){
	 for(i=0;i<N;i++){w[i]=1;}
	 tirages=tiragedetrajectoires(alpha,beta,rho,sigma0,spot,r,q,N,deltat,horizon);
	 G=donnees(tirages,spot,sigma0,r,N,deltat,horizon,T,P,K);
	 poids(w,G,proba,tirages,spot,r,N);
   } else {
	 for(i=0;i<N;i++){w[i]=0;}
	 tirages=tiragedetrajectoires(alpha,beta,rho,sigma0,spot,r,q,N,deltat,horizon);
	 G=donnees(tirages,spot,sigma0,r,N,deltat,horizon,T,P,K);
	 poids(w,G,proba,tirages,spot,r,N);
	 
   }
  
   

  
		  for(j=0;j<N;j++){
			 W=0;
			 for(i=0;i<nu;i++){
			  
			W=W+proba[i]*G[j+N*i];
	
	  
	}
		if(fabs(W)>1E-3){
		  printf("%d %g\n",j,W);
		}
		  }

  free(proba);
  free(G);
  free(tirages);
  free(P);
  free(T);
  free(K);
  free(w);
  
	  return 0;
}
