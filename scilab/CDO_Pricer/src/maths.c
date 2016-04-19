#include        <stdlib.h>
#include        "maths.h"
#include        "math.h"
/** \f$(x-y)_+\f$
 */


/**
 * maximum between two doubles
 *
 * @param x : double
 * @param y : double
 * @return  max(x,y)
 */

double       pp(const double     x,
                       const double     y)
{
  return ( (x > y) ? (x - y) : 0. );
}

/**
 * valuation of a polynom P in a point (be careful : the
 * polynom has non constant term P(0)=0)
 *
 * @param degree : degree of the polynom
 * @param a : array containing the coefficients of the
 * polynom (be careful : a[0] corresponds to the coefficient
 * in front of x^1) 
 * @param x : scalar
 * @return  P(x)
 */

 double       evaluate_poly(const int     degree,
                                  const double  *a,
                                  const double  x)
{
  double      result;
  int         n;

  if (degree == 0) result = 0.;
  else {
    result = a[degree-1] * x;
    for (n = degree-2; n >= 0; n--) 
      result = (result + a[n]) * x;
  }
    
  return (result);
}

/**
 * valuation of the derivative of a polynom P in a point (be careful : the
 * polynom has non constant term P(0)=0)
 *
 * @param degree : degree of the polynom
 * @param a : array containing the coefficients of the
 * polynom (be careful : a[0] corresponds to the coefficient
 * in front of x^1) 
 * @param x : scalar
 * @return  (P')(x)
 */

double       evaluate_dpoly(const int    degree,
                                   const double *a,
                                   const double x)
{
  double      result;
  int         n;

  if (degree == 0) result = 0.;
  else {
    result = degree * a[degree-1];
    for (n = degree-2; n >= 0; n--) 
      result = (n+1.)*a[n] + result * x;
  }
  return (result);
}



static const double     a[] = { -3.969683028665376e+01, 2.209460984245205e+02,
                                -2.759285104469687e+02, 1.383577518672690e+02,
                                -3.066479806614716e+01, 2.506628277459239e+00 };
static const double     b[] = { -5.447609879822406e+01, 1.615858368580409e+02,
                                -1.556989798598866e+02, 6.680131188771972e+01,
                                -1.328068155288572e+01 };
static const double     c[] = { -7.784894002430293e-03, -3.223964580411365e-01,
                                -2.400758277161838e+00, -2.549732539343734e+00, 
                                4.374664141464968e+00, 2.938163982698783e+00 };
static const double     d[] = {  7.784695709041462e-03, 3.224671290700398e-01, 
                                 2.445134137142996e+00, 3.754408661907416e+00 };

/**
 * cumulative normal distribution 
 *
 * @param z : scalar
 * @return N(z) 
 */



double gaussian_cdf(const double z){
if (z>6)    return 0.9999;
 if (z<-6)   return 0.0001;

 double b1=0.31938153;
 double b2=-0.356563782;
 double b3=1.781477937;
 double b4=-1.821255978;
 double b5=1.330274429;
 double c2=0.3989423;
 double p=0.2316419;
 double a=abso(z);
 double t=1/(1+a*p);
 double b=c2*exp((-z)*(z/2));
 double n=((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
 n=1-b*n;
 if (z<0) n=1-n;
 if(z>1) n=0.999;

 return n;
}

/*
double          gaussian_cdf(const double       x)
{
  return ( 1.0 - 0.5 * erfc(x * C_SQRT1_2) );
}

*/

/**
 * inverse of the cumulative normal distribution 
 *
 * @param p : scalar
 * @return  N^{-1}(p)
 */

double          gaussian_inv_cdf(const double       p)
{
  double      q;
  double      r;

  if(p<=0) return -6;

  if ((p < 0.02425)&&(p > 0.0))
    {
      q = sqrt(-2*log(p));
      return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
      ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
 
   if(p>=1) return 6;
 
  else if (p > 0.97575)
    {
      q  = sqrt(-2*log(1-p));
      return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
      ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
  else
    {
      q = p - 0.5;
      r = q*q;
      return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
      (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    }
}

/*

double          gaussian_cdf(const double       x)
{
  return ( 1.0 - 0.5 * erfc(x * C_SQRT1_2) );
}

 Computes the inverse of the cumulative distribution function.
 

double          gaussian_inv_cdf(const double       p)
{
  double      q;
  double      r;

  if (p < 0.02425)
    {
      q = sqrt(-2*log(p));
      return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
      ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
  else if (p > 0.97575)
    {
      q  = sqrt(-2*log(1-p));
      return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
      ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
  else
    {
      q = p - 0.5;
      r = q*q;
      return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
      (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    }
}

*/


double max(double a,double b){
if (a>b) return a;
else     return b;
}




/**LU décomposition **/
 
void ludcmp(double **a,double *index,double d,int n){

/** On doit connaitre à priori le nombre de ligne de la matrice**/

double TINY=1.0e-20;
int i,imax,j,k;
double dp,big,dum,sum,temp;
double *vv;
vv=malloc(n *sizeof(double));
d=1.0;   
for(i=0;i<n;i++){
   big=0.0;
   
   for(j=0;j<n;j++)
    if((temp=abso(a[i][j]))>big) big=temp;
   
    vv[i]=1.0/big;
     
}
  for(j=0;j<n;j++){
    for(i=0;i<j;i++){
     sum=a[i][j];
     
     for(k=0;k<i;k++)
      sum-=a[i][k]*a[k][j];
     
      a[i][j]=sum;
   }
   big=0.0;
  for(i=j;i<n;i++){
    sum=a[i][j];
    for(k=0;k<j;k++) sum-=a[i][k]*a[k][j];
      a[i][j]=sum;
      if((dum=vv[i]*abso(sum))>=big){
      big=dum;
      imax=i;
      }
  }
  if(j!=imax){
    for(k=0;k<n;k++){
     dum=a[imax][k];
     a[imax][k]=a[j][k];
     a[j][k]=dum;
    }
    d=-d;
    vv[imax]=vv[j];
 }
  index[j]=imax;
  if(a[j][j]==0.0) a[j][j]=TINY;
  if(j!=n-1){
  dum=1.0/(a[j][j]);
  for(i=j+1;i<n;i++) a[i][j]*=dum;
  }
 }
}
void lubsk(double **a ,double *index,double *b,int n){
  int i,ii=0,ip,j;
  double sum;
  for(i=0;i<n;i++){
   ip=index[i];
   sum=b[ip];
   b[ip]=b[i];
   if(ii!=0)
    for(j=ii-1;j<i;j++) sum-=a[i][j]*b[j];
   
   else if (sum!=0.0)
   ii=i+1;
   b[i]=sum;
 }
 for(i=n-1;i>=0;i--){
 sum=b[i];
 for(j=i+1;j<n;j++){
   sum-=a[i][j]*b[j];
 }

 b[i]=sum/a[i][i];
 }
}

/**fonction renvoyant l'inverse d'une matrice a n*n  **/

double **inv_mat(double **a,int n){
double *index;
double **y;
int j,i,k;
double *col;

y=malloc(n *sizeof(double*));
for(k=0;k<n;k++){
y[k]=malloc(n *sizeof(double));
}
index=malloc(n*sizeof(double));
col=malloc(n*sizeof(double));
double d;

ludcmp(a,index,d,n);
for(j=0;j<n;j++){
  for(i=0;i<n;i++){
   col[i]=0.0;
  }
 col[j]=1.0;


lubsk(a,index,col,n);
  for(i=0;i<n;i++) y[i][j]=col[i];
}

return y;

}


/** Bessel function I1
 */
static const double     aI1[] = { 0.87890594, 0.51498869, 
                                  0.15084934, 0.02658733, 
                                  0.00301532, 0.00032411 };
static const double     bI1[] = { -0.03988024, -0.00362018,
                                  0.00163801, -0.01031555,
                                  0.02282967, -0.02895312,
                                  0.01787654, -0.00420059 };
double          bessel_I1(const double      x)
{
  double      t = x / 3.75;
  double      it = 3.75 / x;
  double      result;

  t *= t;
  if (x < 3.75) {
    result = x * (0.5 + evaluate_poly(6, aI1, t));
  }
  else {
    result = exp(x) / sqrt(x) * (0.39894228 + evaluate_poly(7, bI1, it));
  }

  return (result);
}

/** Bessel function K1
 */
static const double     aK1[] = {  0.15443144, -0.67278579,
                                   -0.18156897, -0.01919402,
                                   -0.00110404, -0.00004686 };
static const double     bK1[] = {  0.23498619, -0.03655620,
                                   0.01504268, -0.00780353,
                                   0.00325614, -0.00068245 };
double          bessel_K1(const double      x)
{
  double      t = x / 2.;
  double      it = 2. / x;
  double      result;

  t *= t;
  if (x < 2) {
    result = log(x / 2) * bessel_I1(x) + (1 + evaluate_poly(6, aK1,t)) / x;
  }
  else {
    result = exp(-x) / sqrt(x) * (1.25331414 + evaluate_poly(6, bK1, it));
  }

  return (result);
}

/** Generates a normally distributed random number using the polar form of the Box-Muller 
 * transformation and a reject method (very fast and robust). 
 * The global variables \c gaussian_flag, \c first_gaussian and \c second_gaussian are 
 * used.
 */
double          gaussian()
{
  double      x1;
  double      x2;
  double      w;
  static int      gaussian_flag = 0;
  static double   first_gaussian;
  static double   second_gaussian;

  gaussian_flag = (gaussian_flag == 0) ? 1 : 0;
  if (gaussian_flag == 0) return second_gaussian; 
  do {
    x1 = 2.0 * genrand_real3() - 1.0;
    x2 = 2.0 * genrand_real3() - 1.0;
    w = x1 * x1 + x2 * x2;
  } while ( w >= 1.0 );
  w = sqrt( (-2.0 * log(w)) / w );
  first_gaussian = x1 * w;
  second_gaussian = x2 * w; 
    
  return (first_gaussian);
}

double gammas(const double xx){

  double x,y,tmp,ser;
  
static double cof[6]={76.18009172947146,-86.50532032941677,24.014098240824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
int j;
y=x=xx;
tmp=x+5.5;
tmp-=(x+0.5)*log(tmp);
ser=1.000000000190015;
for(j=0;j<=5;j++) {
  y=y+1;
  
ser+=cof[j]/y;
}
return exp(-tmp+log(2.5066282746310005*ser/x));
}


double gammln(double x){
return log(gammas(x));
}

double betacf(double a,double b,double x){
int Maxit=100;
double eps=3.0e-7;
double FPmin=1.0e-30;
int m,m2;
double aa,c,d,del,h,qab,qam,qap;

qab=a+b;
qap=a+1.0;
qam=a-1.0;
c=1.0;
d=1.0-qab*x/qap;
if(abso(d)<FPmin) d=FPmin;
d=1.0/d;
h=d;

for(m=1;m<Maxit;m++){
 m2=2*m;
aa=m*(b-m)*x/((qam+m2)*(a+m2));
d=1.0+aa*d;
if(abso(d)<FPmin) d=FPmin;
c=1.0+aa/c;
if(abso(c)<FPmin) c=FPmin;
d=1.0/d;
h*=d*c;
aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
d=1.0+aa*d;
if(abso(d)<FPmin) d=FPmin;
c=1.0+aa/c;
if(abso(c)<FPmin) c=FPmin;
d=1.0/d;
del=d*c;
h*=del;
if(abso(del-1.0)<=eps) break;
}
return h;
}

double betai(double a,double b,double x){
double bt=0;
if((x<0.0)||(x>1.0)) return 0;
if (x==0||x==1) bt=0;

else 
 bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1-x));

if(x<(a+1.0)/(a+b+2.0)) return bt*(betacf(a,b,x))/a;
else
return 1-bt*(betacf(b,a,1-x))/b;

}

double student_cdf(const double t,const double x){
double pi=3.14159265;
double y, a;                                               
 a=1.0-betai(t/2,1./2,(t/(t+x*x)));

if (x>0)  y=(a*sqrt(pi))/(2*gammas(0.5))+0.5;
else      y=0.5-(a*sqrt(pi))/(2*gammas(0.5));

if(y>0.9999) y=0.9999;

return y;
}


 double  density(const  double t1,const  double x){

return gammas((t1+1)*0.5)/(gammas(t1*0.5)*sqrt(3.14*t1)*exp(((t1+1)/2)*log(1+x*x/t1)));
}

double student_inv_cdf(const double t1,const double x){
int   n=1000;
int i,u,v;
double *p;
double h;
int a;
p=malloc(n*sizeof(double));
h=12./(n-1);
  for(i=0;i<n;i++){
  p[i]=student_cdf(t1,-6+i*h);
  }
 if(x>p[n-1]) return 6;
 if(x<p[0])   return -6;
u=0;
v=n;
i=1;

while((abso(v-u)>1) &&(i<n)){

a=(v+u)/2;
 if (p[a]>x) v=a;
 else u=a;
 

i++;

}
return -6+a*h +((x-p[a])*h)/(p[a+1]-p[a]);
}

double simulate_student(const double t){
double s=0;
int j;
double u;
for(j=0;j<t;j++){
u=gaussian();
s+=u*u;
}
return gaussian()/(sqrt(s/t));
}



 double gauss_density(const double x){
 return (C_SQRT1_2PI  * exp(- 0.5 * x*x) );
}



 double density_chi2(const double t1,const double x){
double a=1/(exp((t1*0.5)*log(2))*gammas(t1*0.5));
if (x>0){
 return a*exp((t1*0.5-1)*log(x))*exp(-x*0.5);
}
else return 0;

}

double abso(const double x){
if(x>0) return x;
else return -x;
}

