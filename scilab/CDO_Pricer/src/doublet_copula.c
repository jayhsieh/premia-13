#include <stdlib.h>
#define _GNU_SOURCE
#include <math.h>
#include "copulas.h"

typedef struct {
    double		rho;
    double		g_rho; 		// g_rho = sqrt(1-rho^2)
    double		u_rho;		// u_rho = rho / g_rho
    double		t1;		//degres de liberte de la premiere student
    double		t2;		//degres de liberte de la deuxieme student
//  double		**cdf;		//la matrice (n+1)*2 pour les points de la cdf et pour calculer son inverse
} doublet_params;


static double		**cdf;		//la matrice (n+1)*2 pour les points de la cdf et pour calculer son inverse
static double		t1_student_density;
static double		t2_student_density;


static const double 	dp[] =		//parametres pour la cdf de la doublet (bornes integration)
{ -3e+1,
  3e+1,
  -1e+2,
  1e+2
};

static const int	di[] = 		//parametres pour l'integration - nombres de points
{ 400,
  400
};



static double		student_density(const copula		*cop,
					const double		x)
{
    double t;
    t=t1_student_density;  
    return ( tgamma((t+1)/2)/(tgamma(t/2)*sqrt(M_PI * t)*pow(1+x*x/t,(t+1)/2)) ); //fdr de la loi student a t deg de lib
}


double 			betai(double a,
			      double b,
			      double x) 		//Returns the incomplete beta function Ix(a, b).
{
double betacf(double a, double b, double x);
double gammln(double xx);
double bt;
if (x == 0.0 || x == 1.0) bt=0.0;
else 					//Factors in front of the continued fraction.
bt=exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1.0-x));
if (x < (a+1.0)/(a+b+2.0)) 		//Use continued fraction directly.
return bt*betacf(a,b,x)/a;
else 					//Use continued fraction after making the symmetry transformation
return (1.0-bt*betacf(b,a,1.0-x)/b); 
}


double			 betacf(double a,
				double b,
				double x)

{
int m,m2;
double aa,c,d,del,h,qab,qam,qap,MAXIT,EPS,FPMIN;
MAXIT=100;
EPS=3.0e-7;
FPMIN=1.0e-30;
qab=a+b;
qap=a+1.0;
qam=a-1.0;
c=1.0;
d=1.0-qab*x/qap;
if (fabs(d) < FPMIN) d=FPMIN;
d=1.0/d;
h=d;
for (m=1;m<=MAXIT;m++) {
m2=2*m;
aa=m*(b-m)*x/((qam+m2)*(a+m2));
d=1.0+aa*d;
if (fabs(d) < FPMIN) d=FPMIN;
c=1.0+aa/c;
if (fabs(c) < FPMIN) c=FPMIN;
d=1.0/d;
h *= d*c;
aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
d=1.0+aa*d;
if (fabs(d) < FPMIN) d=FPMIN;
c=1.0+aa/c;
if (fabs(c) < FPMIN) c=FPMIN;
d=1.0/d;
del=d*c;
h *= del;
if (fabs(del-1.0) < EPS) break;
}
return h;
}


static double		student_cdf(const double		x)
{
   double t;
   double y;
   t=t2_student_density;
   if(x>=0) y=1-1/2*betai(t/2,1/2,t/(t+x*x));
   else y=1/2*betai(t/2,1/2,t/(t+x*x));
   return(y);
}

/* Inverse doublet cumulative distribution function */

static double		doublet_inv_cdf(const double 		x)
{
    int 		n=di[0];
    int			i=2;
    int			d=2;
    int			u1=0;			//borne inf de la dichotomie
    int			u2=(int) (floor(n/2));	//borne sup de la dichotomie
    
    while(i<n & abs(u2-u1)>1){
    	if(cdf[1][u2]>x) u2-=(int)(floor((n+1)*pow(1/2,i)));
	else u1+=(int)(floor((n+1)*pow(1/2,i)));
	i++;
    }
    
    return(cdf[0][u2]+(x-cdf[1][u2])*(cdf[0][u2+1]-cdf[0][u2])/(cdf[1][u2+1]-cdf[1][u2]));
}

static double		*doublet_compute_prob(const copula		*gc,
					       const double		f_t,
					       const grid		*v)
{
    double		*result;
    doublet_params	*p;
    double		a;
    int	i;
    
    p = gc->parameters;
    result = malloc(cop->size * sizeof(double));
    a = sqrt(p->t2 / (p->t2 - 2)) * doublet_inv_cdf(f_t) / p->g_rho;
    for (i = 0; i < cop->size; i++) {
        result[i] = student_cdf(a - sqrt(p->t2 / (p->t2 - 2) * (p->t1 - 2) / p->t1) *
p->u_rho * cop->points[i]);
//      printf("%g\t%g\n", cop->points[i], result[i]);
    }
    
    return (result);
}


//fonction pour la double integration de la cdf
static double		f_cdf(const double rho,
			      const double t1,
			      const double t2,
			      double x1,
			      double x2)
{
    double u;
    double v;
    u = pow(1 + (x1 - x2)*(x1 - x2) / (rho*rho*(t1-2)) , (t1+1)/2);
    v = pow(1 + x2*x2 / ((1-rho*rho)*(t2-2)) , (t2+1)/2);
    return (1./(u*v));
}

void			init_cdf_doublet(const double rho,
					   const double t1,
					   const double t2,
					   int n,
					   int m)
{
//  double 		**cdf_doublet;
    double		coef1;			//coef pour la cdf
    double		coef2;			//coef pour la deuxieme integrale
    int 		i=0;
    int			j=0;
    double		a1;			//borne inf de la cdf
    double		b1;			//borne sup de la cdf
    double		a2;			//borne inf de la deuxieme integrale
    double		b2;			//borne sup de la deuxieme integrale
    double		s1;			//variable de stockage
    double		s2;			//autre variable de stockage
    
    
    cdf=malloc(2*sizeof(double*));		//on declare le tableau pour la cdf
    cdf[0]=malloc((n+1)*sizeof(double));	//on declare la colonne des points
    cdf[1]=malloc((n+1)*sizeof(double));	//on declare la colonne des valeurs aux points
    coef1 = (tgamma((t1+1)/2) * tgamma((t2+1)/2));
    coef1 = coef1 / (tgamma(t1/2) * tgamma(t2/2) * M_PI * rho * sqrt((1 - rho*rho) * (t1 - 2) * (t2 - 2))); 
    a1=dp[0];
    b1=dp[1];
    a2=dp[2];
    b2=dp[3];
    coef1 = coef1*(b1 - a1)/n;			//methode des trapezes
    coef2 = (b2 - a2)/m;				//methode des trapezes
    printf("test : f_cdf = %f \n",pow(a1,2));
    s2=0;    
    printf("avant boucle\n");
    for(i=0;i<n+1;i++){
    cdf[0][i]= a1 + (b1 - a1)*i/( (double )n); 		//point considere = x1;
    s1=0;
    s1+=f_cdf(rho,t1,t2,cdf[0][i],a2)/2;
    s1+=f_cdf(rho,t1,t2,cdf[0][i],b2)/2;
    	for(j=1;j<m;j++)
	s1+=f_cdf(rho,t1,t2,cdf[0][i],a2 + (b2 - a2)*j/((double) m));	//methode des trapezes
    if (i>0) {
  //  printf("s1 = %f, s2= %f ,\n",s1,s2);
    cdf[1][i]=(s1+s2)/2*coef1*coef2+cdf[1][i-1];		//methode des trapezes
  //  printf("X=%g, Y=%g, i=%f\n", cdf[1][i], cdf[1][i-1], cdf[0][i]);
    }
    else
    cdf[1][i]=0;						//premiere boucle, initialisation
    
    s2=s1;
    }	
  //  printf("ok %f, %f, %f\n",cdf[0][n],cdf[1][n],cdf[1][50]);
    return;
}

copula			*init_doublet_copula(const double	rho,
						 const double	t1,
						 const double	t2)
{
    copula		*gc;
    doublet_params	*p;
 // double		**cdf;			//on n'initialise pas comme c'est defa fait en statique
    int			n=di[0];
    int			m=di[1];
    
    
    gc = malloc(sizeof(copula));
    gc->name = "one factor double-t copula";
    t1_student_density=t1; 			//on se sert de t1 car c'est celui qui va servir pour l'integration
    t2_student_density=t2;
    gc->density = student_density; 			//la densite est celle d'une student a t1 degres de liberte
//  cdf=malloc(2*sizeof(double*));	//on declare le tableau pour la cdf
//  cdf[0]=malloc((n+1)*sizeof(double));	//on declare la colonne des points
//  cdf[1]=malloc((n+1)*sizeof(double));	//on declare la colonne des valeurs aux points
    printf("avant init\n");
    init_cdf_doublet(rho, t1, t2, n, m);	//initialise le tableau **cdf
    printf("avant init\n");
//    gc->phi = student_cdf;			//la cdf 
//    gc->inv_phi = doublet_inv_cdf;
    gc->compute_cond_prob = doublet_compute_prob;
    p = malloc(sizeof(doublet_params));
    p->rho = rho;
    p->g_rho = sqrt(1.0 - rho*rho);
    p->u_rho = rho / p->g_rho;
    p->t1 = t1;
    p->t2 = t2;
    gc->parameters = p;
    printf("init ok\n");

    return (gc);
}


