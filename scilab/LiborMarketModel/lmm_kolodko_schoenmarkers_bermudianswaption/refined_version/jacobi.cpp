#include <math.h>

#include <iostream>
#include <valarray>
#include "nrutil.h"
#include "nrutil.c"

using namespace std;

class Matrice
{
 public:
  Matrice(int _n, double _phi, int _d);
  Matrice();

  int n; //dimnesion of matrix ro
  int d; //rank
  double phi;
  valarray<valarray<double> > a; //matrix ro 
  valarray<double> spec; //ro spectrum
  valarray<valarray<double> > v; // eigenvectors
  valarray<valarray<double> > les_ei; //ei vectors
  valarray<valarray<double> > ro_d; //reduced matrix
  
  void jacobi();//cf Numerical recipe
  void eigsrt();//cf Numerical recipe
  void reduce();
};

#define ROTATE(a,i,j,k,l) g=a[i-1][j-1];h=a[k-1][l-1];a[i-1][j-1]=g-s*(h+g*tau); a[k-1][l-1]=h+s*(g-h*tau);

Matrice::Matrice(int _n, double _phi,int _d):n(_n),d(_d),phi(_phi)
{  
  spec.resize(n); 
  int iq,ip;
  v.resize(n); 
  for (ip=1;ip<=n;ip++) { //Initialize to the identity matrix.
    v[ip-1].resize(n);
    for (iq=1;iq<=n;iq++) v[ip-1][iq-1]=0.0;
    v[ip-1][ip-1]=1.0;
  }  
  a.resize(n);
  for (ip=1;ip<=n;ip++) 
    {
      a[ip-1].resize(n);  
      for (iq=1;iq<=n;iq++) 
	{
	  a[ip-1][iq-1] = exp(-phi*abs(ip-iq));
	}
    }
  les_ei.resize(n); 
  for (ip=1;ip<=n;ip++) { //Initialize to the identity matrix.
    les_ei[ip-1].resize(d);
  }  
  ro_d.resize(n); 
  for (ip=1;ip<=n;ip++) { //Initialize to the identity matrix.
    ro_d[ip-1].resize(n);
  }  
}

Matrice::Matrice()//nada
{}

void Matrice::jacobi()
/* Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n]. On
   output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
   v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
   a. nrot returns the number of Jacobi rotations that were required. */
{
  int *nrot;
  int j,iq,ip,i;
  float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
  b=vector(1,n);
  z=vector(1,n);
  for (ip=1;ip<=n;ip++) { //Initialize b and d to the diagonal
    b[ip]=spec[ip-1]=a[ip-1][ip-1]; //of a.
    z[ip]=0.0; //This vector will accumulate terms of the form tapq as in equation (11.1.14).
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++) { //Sum off-diagonal elements.
      for (iq=ip+1;iq<=n;iq++)
	sm += fabs(a[ip-1][iq-1]);
    }
    if (sm == 0.0) { //The normal return, which relies on quadratic convergence to machine underflow.
      free_vector(z,1,n);
      free_vector(b,1,n);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n); //...on the first three sweeps.
    else
      tresh=0.0; //...thereafter.
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
	g=100.0*fabs(a[ip-1][iq-1]); //After four sweeps, skip the rotation if the off-diagonal element is small.
	if (i > 4 && (float)(fabs(spec[ip-1])+g) == (float)fabs(spec[ip-1]) && (float)(fabs(spec[iq-1])+g) == (float)fabs(spec[iq-1])) a[ip-1][iq-1]=0.0;
	else if (fabs(a[ip-1][iq-1]) > tresh) {
	  h=spec[iq-1]-spec[ip-1];
	  if ((float)(fabs(h)+g) == (float)fabs(h)) t=(a[ip-1][iq-1])/h; //t = 1/theta
	  else {
	    theta=0.5*h/(a[ip-1][iq-1]); //Equation (11.1.10).
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip-1][iq-1];
	  z[ip] -= h;
	  z[iq] += h;
	  spec[ip-1] -= h;
	  spec[iq-1] += h;
	  a[ip-1][iq-1]=0.0;
	  for (j=1;j<=ip-1;j++) { //Case of rotations 1 \u2264 j < p.
	    ROTATE(a,j,ip,j,iq)
	      }
	  for (j=ip+1;j<=iq-1;j++) { //Case of rotations p < j < q.
	    ROTATE(a,ip,j,j,iq)
	      }
	  for (j=iq+1;j<=n;j++) { //Case of rotations q < j \u2264 n.
	    ROTATE(a,ip,j,iq,j)
	      }
	  for (j=1;j<=n;j++) {
	    ROTATE(v,j,ip,j,iq)
	      }
	  ++(*nrot);
	}
      }
    }
    for (ip=1;ip<=n;ip++) {
      b[ip] += z[ip];
      spec[ip-1]=b[ip]; //Update d with the sum of tapq,
      z[ip]=0.0; //and reinitialize z.
    }
  }
  nrerror("Too many iterations in routine jacobi");
}

void Matrice::eigsrt()
/*
  Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output from jacobi
  (§11.1) or tqli (§11.3), this routine sorts the eigenvalues into descending order, and rearranges
  the columns of v correspondingly. The method is straight insertion.*/
{
  int k,j,i;
  float p;
  for (i=0;i<n-1;i++) {
    p=spec[k=i];
    for (j=i;j<n;j++)
      if (spec[j] >= p) p=spec[k=j];
    if (k != i) {
      spec[k]=spec[i];
      spec[i]=p;
      for (j=0;j<n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
}	

void Matrice::reduce()
{   
  cout << endl;
  jacobi ();
  eigsrt ();
  int ip, iq, j;
  double norme;
  for (ip=1;ip<=n;ip++) 
    {  
      for (iq=1;iq<=d;iq++) 
	{
	  v[ip-1][iq-1] *= sqrt(spec[iq-1]); 
	}
    }  
  for (ip=1;ip<=n;ip++) 
    {  
      norme =0.;
      for (iq=1;iq<=d;iq++) 
	{
	  norme += v[ip-1][iq-1]*v[ip-1][iq-1];
	}
      norme = sqrt(norme);
      for (iq=1;iq<=d;iq++) 
	{
	  les_ei[ip-1][iq-1] = v[ip-1][iq-1]/norme;
	}
    } 
  for (ip=1;ip<=n;ip++) 
    {  
      for (iq=1;iq<=n;iq++) 
	{
	  ro_d[ip-1][iq-1] = 0.;
	  for (j=1;j<=d;j++) 
	    {
	      ro_d[ip-1][iq-1] += les_ei[ip-1][j-1]*les_ei[iq-1][j-1];
	    }
	}
    }
}
