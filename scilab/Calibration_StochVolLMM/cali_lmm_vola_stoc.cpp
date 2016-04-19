#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <vector>
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <complex>
#include <stddef.h>
#include <stdlib.h>

#include "intg.h"
#include "realfft.h"
#include "normal_df.h"
#include "min.h"

using namespace std;
using std::string;
#define FREE_ARG char*

double * T, * zero_bond, * forward_rate, * lf,  * ke, * sigma, * rho,* strike, ** caplet_vola, ** caplet_price1,  ** gamma1, ** gamma2, r,C1, C2, A, B, G, * ifftresult , * ifftresult1, * ifftresult2, * ifftresult3, * ifftresult4, * ifftresult5, * xi, * vxi, * vxi1, *vxi2, * vxi3, * vxi4, * vxi5, r0 ;
int TN, SN, maturity;
int kmax=8*512;//the step number of inverse fourier transform //
double h=0.03;//the step size of inverse fourier transform
FILE  *output=fopen ("cali_result", "w"); // result output file


double *allocate_vector(int nl, int nh)//allocate a double vector with subscript range v[nl..nh]
{
	double *v;
	v=(double *)malloc((size_t) ((nh-nl+1+1)*sizeof(double)));
	if (!v)		printf("allocation failure in vector()\n");
	return v-nl+1;
}
void free_vector(double *v, int nl, int nh)//free a double vector allocated with allocate_vector()
{
	free((FREE_ARG)(v+nl-1));
}
double ** allocate_matrix(int nrl, int nrh, int ncl, int nch) //allocate a double matrix with subscript range v[nl..nh]
{
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	double **m;
	
	m=(double **) malloc ((size_t)((nrow+1)*sizeof(double*)));
	if (!m)		printf("allocation failure 1 in matrix()\n");
	m+=1;
	m-=nrl;

	m[nrl]=(double *) malloc ((size_t)((nrow*ncol+1)*sizeof(double)));
	if (!m[nrl])	printf("allocation failure 2 in matrix()\n");
	m[nrl] += 1;
	m[nrl] -= ncl;

	for (i=nrl+1; i<= nrh; i++) m[i]=m[i-1]+ncol;

	return m;
}
void free_matrix(double **m, int nrl, int nrh, int ncl, int nch) //free a double matrix with subscript range v[nl..nh]
{
	free((FREE_ARG) (m[nrl]+ncl-1));
	free((FREE_ARG) (m+nrl-1));
}
void allocate_memory(int tenor_number, int strike_number)
{
	T=allocate_vector(0, tenor_number);
	zero_bond=allocate_vector(0, tenor_number);
	forward_rate=allocate_vector(0, tenor_number);
	lf=allocate_vector(0, tenor_number);
	ke=allocate_vector(0, tenor_number);
	sigma=allocate_vector(0, tenor_number);
	rho=allocate_vector(0, tenor_number);
	strike=allocate_vector(0, strike_number);	
	caplet_vola= allocate_matrix(0,strike_number, 0, tenor_number);
	caplet_price1=allocate_matrix(0,strike_number, 0, tenor_number);
	gamma1=allocate_matrix(0, tenor_number,0, tenor_number);
	gamma2=allocate_matrix(0, tenor_number,0, tenor_number);
	ifftresult=allocate_vector(0, 2*kmax+1);
	ifftresult1=allocate_vector(0, 2*kmax+1);
	ifftresult2=allocate_vector(0, 2*kmax+1);
	ifftresult3=allocate_vector(0, 2*kmax+1);
	ifftresult4=allocate_vector(0, 2*kmax+1);
	ifftresult5=allocate_vector(0, 2*kmax+1);
	xi=allocate_vector(0, 2*kmax+1);
	vxi=allocate_vector(0, 2*kmax+1);
	vxi1=allocate_vector(0, 2*kmax+1);
	vxi2=allocate_vector(0, 2*kmax+1);
	vxi3=allocate_vector(0, 2*kmax+1);
	vxi4=allocate_vector(0, 2*kmax+1);
	vxi5=allocate_vector(0, 2*kmax+1);
}
void free_memory(int tenor_number, int strike_number)
{
	free_vector(T, 0, tenor_number);
	free_vector(zero_bond, 0, tenor_number);
	free_vector(forward_rate, 0, tenor_number);
	free_vector(strike, 0, strike_number);
	free_matrix(caplet_vola, 0, strike_number, 0, tenor_number);
	free_matrix(caplet_price1,  0, strike_number, 0, tenor_number);
	free_vector(lf, 0, tenor_number);
	free_vector(ke, 0, tenor_number);
	free_vector(sigma, 0, tenor_number);
	free_vector(rho, 0, tenor_number);
	free_matrix(gamma1, 0, tenor_number, 0, tenor_number);
	free_matrix(gamma2,  0, tenor_number, 0, tenor_number);
	free_vector(ifftresult, 0, 2*kmax+1 );
	free_vector(ifftresult1, 0, 2*kmax+1 );
	free_vector(ifftresult2, 0, 2*kmax+1 );
	free_vector(ifftresult3, 0, 2*kmax+1 );
	free_vector(ifftresult4, 0, 2*kmax+1 );
	free_vector(ifftresult5, 0, 2*kmax+1 );
	free_vector(xi, 0, 2*kmax+1 );
	free_vector(vxi, 0, 2*kmax+1 );
	free_vector(vxi1, 0, 2*kmax+1 );
	free_vector(vxi2, 0, 2*kmax+1 );
	free_vector(vxi3, 0, 2*kmax+1 );
	free_vector(vxi4, 0, 2*kmax+1 );
	free_vector(vxi5, 0, 2*kmax+1 );
} 
void initialize(int tenor_number, int strike_number)
{
	for(int i=0; i< tenor_number; i++)
	{
		T[i]=0.0;
		zero_bond[i]=0.0;
		forward_rate[i]=0.0;
		lf[i]=0.0;
		ke[i]=0.0;
		sigma[i]=0.0;
		rho[i]=0.0;
		for (int k=0; k<strike_number; k++)
		{
			strike[k]=0.0;
			caplet_vola[k][i]=0.0;
			caplet_price1[k][i]=0.0;
		}
		for (int j=0; j<tenor_number; j++)
		{
			gamma1[i][j]=0.0;
			gamma2[i][j]=0.0;
		}
	}
	for (int n=0; n<2*kmax+1; n++)
	{
		ifftresult[n]=0.0;
		ifftresult1[n]=0.0;
		ifftresult2[n]=0.0;
		ifftresult3[n]=0.0;
		ifftresult4[n]=0.0;
		ifftresult5[n]=0.0;
		xi[n]=0.0;
		vxi[n]=0.0;
		vxi1[n]=0.0;
		vxi2[n]=0.0;
		vxi3[n]=0.0;
		vxi4[n]=0.0;
		vxi5[n]=0.0;
	}
}

void InitialCurveAndZeroBond(int tenor_number, double tenor)
{
	T[0]=0.0;
	zero_bond[0]=1.0;
	FILE * input = fopen ("initial_curve.dat", "r");
	int tenor_order;
	float a1;        //temperory variance for saving datas from file	
	if (input !=NULL)
	{
		while (! feof (input) )
 		{
			fscanf (input, "%f", &a1);
			tenor_order=int(a1/tenor);
			if (tenor_order <=tenor_number && (a1-tenor_order*tenor)==0.0 )
			{
					fscanf (input, "%f", &a1);
					forward_rate[tenor_order]=a1;
					T[tenor_order]=tenor_order*tenor;
					zero_bond[tenor_order]=zero_bond[tenor_order-1]/(1.+(T[tenor_order]-T[tenor_order-1])*forward_rate[tenor_order]);
  			}
			else
			{
 				break;
			}
 		}
 		fclose (input);
		T[tenor_number-1]=(tenor_number-1)*tenor;
	}
}
void ReadMarketData(int tenor_number, double tenor, int strike_number) //Read data from file "market_data.txt" 
{
	float a1;        //temperory variance for saving datas from file
	int tenor_order; 
	FILE *input = fopen ("cap_strike_matrix.dat","r");
	for (int i=0; i<tenor_number; i++)
	{
		for(int k=0; k< strike_number;k++)
		{
			caplet_vola[k][i]=0.0;
		}
	}
	if (input!= NULL)
	{
		for(int i=1; i< strike_number; i++)
		{
			fscanf(input, "%f", &a1 );
			strike[i]=a1;
		}
		while (! feof (input) )
 		{
			fscanf (input, "%f", &a1);
			tenor_order=int(a1/tenor);
			if (tenor_order <=tenor_number && (a1-tenor_order*tenor)==0.0 )
			{
				for(int i=1; i<strike_number; i++)
				{
					fscanf (input, "%f", &a1);
					caplet_vola[i][tenor_order]=a1;
				}
			}
			else
			{
				break;
			}
		}
		fclose (input);
	}
	else cout << "Unable to open file";
}


double Bl(double k, double x, double vola, double time) // black formula 
{
	double d1=(log(x/k)+vola*vola*time/2.)/(vola*sqrt(time));
	double d2=d1-vola*sqrt(time);
	return double(x*normal_df(d1)-k*normal_df(d2));
}
void CapletPrice1(int tenor_number, int strike_number) //caplet price from market volas using black formula
{
	for(int i=1; i< tenor_number; i++)
	{
		for(int j=1; j< strike_number; j++)
		{
			if(caplet_vola[j][i]!=0.0)
			{
 				caplet_price1[j][i]=zero_bond[i]*(T[i+1]-T[i])*Bl(strike[j], forward_rate[i], caplet_vola[j][i], T[i]);
			}
			else caplet_price1[j][i]=0.0;
		}
	}
}
double g_intg(double a, double b, double g, double T) // integral of g(T-s)*g(T-s) from 0 to T, g(s)=g+(1-g+a*s)*exp(-b*s)
{
	return double (-(-0.4e1 * g * g * pow(b, 0.3e1) * T + 0.8e1 * exp(-b * T) * g * b * b - 0.8e1 * exp(-b * T) * g * g * b * b + 0.8e1 * g * a * exp(-b * T) * b * b * T + 0.8e1 * g * a * exp(-b * T) * b + 0.2e1 * exp(-0.2e1 * b * T) * b * b - 0.4e1 * exp(-0.2e1 * b * T) * g * b * b + 0.4e1 * a * exp(-0.2e1 * b * T) * b * b * T + 0.2e1 * a * exp(-0.2e1 * b * T) * b + 0.2e1 * exp(-0.2e1 * b * T) * g * g * b * b - 0.4e1 * g * a * exp(-0.2e1 * b * T) * b * b * T - 0.2e1 * g * a * exp(-0.2e1 * b * T) * b + 0.2e1 * a * a * exp(-0.2e1 * b * T) * b * b * T * T + 0.2e1 * a * a * exp(-0.2e1 * b * T) * b * T + a * a * exp(-0.2e1 * b * T) - 0.4e1 * g * b * b + 0.6e1 * g * g * b * b - 0.6e1 * g * a * b - 0.2e1 * b * b - 0.2e1 * b * a - a * a) * pow(b, -0.3e1) / 0.4e1);
}
double Gamma2(int i, int j, double c)// as given by (13)of page 7 in Paper of Belomestny, Mathew and Schoenmakers
{
	double a=0.;
	if (j>=i)  a=c*sqrt(g_intg(A, B, G, T[i])/T[i]);
	else a=0.;
	return a;
}


void kappa(double r0, double k0, double sigma0, double rho0, double lf0, int tenor_number, int k, int j,double result[])
{
	double sum=0.0;
	for(int i=j+1; i<tenor_number-1; i++)
	{
		sum=sum+(T[i+1]-T[i])*forward_rate[i]*Gamma2(i,k,lf[i])/(1.+(T[i+1]-T[i])*forward_rate[i]);
	}
	result[0]= k0-r*sigma0*rho0*sum;
	result[1]=-sigma0*rho0*sum;
	result[2]=1.;
	result[3]=-r*rho0*sum;
	result[4]=-r*sigma0*sum;
	result[5]=0.0;
}
void beta (double r,double sigma0, double lf0, int k,int j,double result[])
{	
	result[0]= r*Gamma2(j, k, lf0)*sigma0;
	result[1]=Gamma2(j, k, lf0)*sigma0;
	result[2]=0.0;
	result[3]=r*Gamma2(j,k,lf0);
	result[4]=0.0;
	result[5]=r*sigma0*sqrt(g_intg(A, B, G, T[j])/T[j]);
}
void alpha(double Kappa[] ,double b[],double rho0, double z, complex<double> result[])//d_{j,k} as given by (38) in page 17 of paper
{
	for (int i=0; i<6; i++)
	{
		if(i==4)	result[i]=Kappa[i]-I*b[0]*(z-I);
		else		result[i]=Kappa[i]-I*b[i]*rho0*(z-I);
	}
}
void fd(complex<double> a[] ,double b[], double z, complex<double> result[])//d_{j,k} as given by (38) in page 17 of paper
{
	complex<double> temp=a[0]*a[0]+b[0]*b[0]*(z*z-I*z);
	if(real(temp)*real(temp)+imag(temp)*imag(temp)>1.e-10)
	{	
		result[0]=sqrt(a[0]*a[0]+b[0]*b[0]*(z*z-I*z));
		for (int i=1; i<6; i++)
		{
			result[i]= (a[0]*a[i]+b[0]*b[i]*(z*z-I*z))/result[0];
		}
	}
	else
	{
		for (int i=0; i<6; i++)
		{
			result[i]=0.0;
		}
	}
}
void fg(complex<double> a[] ,complex<double> d[],complex<double> result[])//g_{j,k} as given by (38) in page 17 of paper
{
	if((real(a[0]-d[0])*real(a[0]-d[0])+imag(a[0]-d[0])*imag(a[0]-d[0]))<= 1.0e-10)
	{	
		for (int i=0; i<6; i++)
		{
			result[i]=0.0;
		}
 	}
	else	
	{	
		result[0]=(a[0]+d[0])/(a[0]-d[0]);
		for (int i=1; i<6; i++)
		{	
     			result[i]=((a[i]+d[i])*(a[0]-d[0])-(a[0]+d[0])*(a[i]-d[i]))/(a[0]-d[0])/(a[0]-d[0]);
		}
	}
}
void fA(complex<double> a[] ,double b[],complex<double> d[],complex<double> g[], double k0, double sigma0,  int j, complex<double> result[])//A_{j,k} as given by (38) in page 17 of paper
{
	complex<double> temp1 = pow((1.0-g[0]),2.0);
	complex<double> temp2 = exp(-d[0]*T[j])-g[0];
	if (real(temp1)*real(temp1)+imag(temp1)*imag(temp1)>1.e-10 && real(temp2)*real(temp2)+imag(temp2)*imag(temp2)> 1.e-10)
	{
		result[0]=k0*pow(sigma0,-2.)*((a[0]-d[0])*T[j]-2.0*log((exp(-d[0]*T[j])-g[0])/(1.-g[0])));
		for (int i=1; i< 6; i++)
		{
			if(i==2) 
			{
				result[i]=pow(sigma0,-2.)*((a[0]-d[0])*T[j]-2.0*log((exp(-d[0]*T[j])-g[0])/(1.-g[0])))+k0*pow(sigma0,-2.)*((a[i]-d[i])*T[j]-2.0*(1.-g[0])/(exp(-d[0]*T[j])-g[0])*((-d[i]*T[j]*exp(-d[0]*T[j])-g[i])/(1.-g[0])+(exp(-d[0]*T[j])-g[0])*g[i]/pow((1.-g[0]),2.)));
			}
			else if(i==3)
			{
				result[i]=k0*(-2.*pow(sigma0,-3.)*((a[0]-d[0])*T[j]-2.0*log((exp(-d[0]*T[j])-g[0])/(1.-g[0])))+pow(sigma0,-2.)*((a[i]-d[i])*T[j]-2.0*(1.-g[0])/(exp(-d[0]*T[j])-g[0])*((-d[i]*T[j]*exp(-d[0]*T[j])-g[i])/(1.-g[0])+(exp(-d[0]*T[j])-g[0])*g[i]/pow((1.-g[0]),2.))));
			}
			else
			{
				result[i]=k0*pow(sigma0, -2.)*((a[i]-d[i])*T[j]-2.0*(1.-g[0])/(exp(-d[0]*T[j])-g[0])*((-d[i]*T[j]*exp(-d[0]*T[j])-g[i])/(1.-g[0])+(exp(-d[0]*T[j])-g[0])*g[i]/pow((1.-g[0]),2.)));
			}
		}
  	}
	else
	{
		for(int i=0; i<6; i++)
		{
			result[i]=0.0;
		}
 	}
}
void fB(complex<double> a[] ,double b[],complex<double> d[],complex<double> g[],double sigma0,  int j, complex<double> result[])//A_{j,k} as given by (38) in page 17 of paper
{
	if (pow(real(pow(sigma0, 0.3e1) * pow(exp(-d[0] * T[j]) - g[0], 0.2e1)),2.0)+pow(imag(pow(sigma0, 0.3e1) * pow(exp(-d[0] * T[j]) - g[0], 0.2e1)),2.0) >= 1.0e-10)	
	{
		result[0]=(a[0]+d[0])*(exp(-d[0]*T[j])-1.0)/(sigma0*sigma0*(exp(-d[0]*T[j])-g[0]));
		for (int i=1; i<6; i++)
		{
			if (i==3)	result[i]=(a[i] + d[i]) * (exp(-d[0] * T[j]) - 0.1e1) * pow(sigma0, -0.2e1) / (exp(-d[0] * T[j]) - g[0]) - (a[0] + d[0]) * d[i] * T[j] * exp(-d[0] * T[j]) * pow(sigma0, -0.2e1) / (exp(-d[0] * T[j]) - g[0]) - 0.2e1 * (a[0] + d[0]) * (exp(-d[0] * T[j]) - 0.1e1) * pow(sigma0, -0.3e1) / (exp(-d[0] * T[j]) - g[0]) - (a[0] + d[0]) * (exp(-d[0] * T[j]) - 0.1e1) * pow(sigma0, -0.2e1) * pow(exp(-d[0] * T[j]) - g[0], -0.2e1) * (-d[i] * T[j] * exp(-d[0] * T[j]) - g[i]);
	 		else	result[i]=(a[i] + d[i]) * (exp(-d[0] * T[j]) - 0.1e1) * pow(sigma0, -0.2e1) / (exp(-d[0] * T[j]) - g[0]) - (a[0] + d[0]) * d[i] * T[j] * exp(-d[0] * T[j]) * pow(sigma0, -0.2e1) / (exp(-d[0] * T[j]) - g[0]) - (a[0] + d[0]) * (exp(-d[0] * T[j]) - 0.1e1) * pow(sigma0, -0.2e1) * pow(exp(-d[0] * T[j]) - g[0], -0.2e1) * (-d[i] * T[j] * exp(-d[0] * T[j]) - g[i]);
		}
	}
	else
	{
		for(int i=0; i<6; i++)
		{
			result[i]=0.0;
		}
	}
}
void Phi(double r0, double k0, double sigma0, double rho0, double lf0, int tenor_number, int j, double z, complex<double> result[]) //phi_{j+1}(z-I;T_j,l,v) as given by (38) in page 17 of paper
{
	r=r0;
	lf[j]=lf0;

	double re [6];
	double im [6];
	double Kappa [6];
	double b [6];
	complex<double> a [6];
	complex<double> d [6];
	complex<double> g [6];
	complex<double> fa [6];
	complex<double> fb [6];
	complex<double> test [6];
	complex<double> testde[6];

	re[0]=-0.5 *(1.-r*r)* lf[j]*lf[j]*g_intg(A, B, G, T[j]) * z * z;
	im[0]= 0.5 *(1.-r*r)* lf[j]*lf[j]*g_intg(A, B, G, T[j]) *z;
	re[1]= r*lf[j]*lf[j]*g_intg(A, B, G, T[j]) * z * z;
	im[1]= -r*lf[j]*lf[j]*g_intg(A, B, G, T[j]) *z;
	re[5]=-(1.-r*r)*lf[j]*g_intg(A, B, G, T[j]) * z * z;
	im[5]=(1.-r*r)* lf[j]*g_intg(A, B, G, T[j]) *z;
	for (int i=2; i<5; i++)
	{
		re[i]=0.0;
		im[i]=0.0;
	}

	ke[j]=k0;
	sigma[j]=sigma0;
	rho[j]=rho0;
	kappa(r, ke[j], sigma[j], rho[j], lf[j], tenor_number, j, j,Kappa);
	beta(r, sigma[j],lf[j], j, j, b);
	alpha(Kappa,b,rho[j], z,a);
	fd(a,b,z,d);
	fg(a,d,g);
	fA(a,b,d,g, ke[j], sigma[j],j,fa);
	fB(a ,b, d, g, sigma[j], j, fb);
	for (int i=0; i<6; i++)
	{	
		re[i]=re[i]+real(fa[i])+real(fb[i]);
		im[i]=im[i]+imag(fa[i])+imag(fb[i]);
	}
	for (int k=j+1; k<tenor_number-1; k++)
	{
		kappa(r, ke[k], sigma[k], rho[k], lf[j], tenor_number, k, j,Kappa);
		beta(r, sigma[k],lf[j], k, j, b);
		alpha(Kappa,b,rho[k], z,a);
		fd(a,b,z,d);
		fg(a,d,g);
		fA(a,b,d,g, ke[k], sigma[k],j,fa);
		fB(a ,b, d, g, sigma[k], j, fb);
		for (int i=0; i<6; i++)
		{	
			re[i]=re[i]+real(fa[i])+real(fb[i]);
			im[i]=im[i]+imag(fa[i])+imag(fb[i]);
		}
	}

	result[0]=exp(re[0])*(cos(im[0])+I*sin(im[0]));
	for (int i=1; i<6; i++)
	{
		result[i]=re[i]*exp(re[0])*cos(im[0])-im[i]*exp(re[0])*sin(im[0])
		+I*(re[i]*exp(re[0])*sin(im[0])+im[i]*exp(re[0])*cos(im[0]));
	}

	test[0]=exp(-0.5 *(1.-r*r)* lf[j]*lf[j]*g_intg(A, B, G, T[j]) * (z * z- I*z));
	testde[0]=0.0;
	testde[1]= r* lf[j]*lf[j]*g_intg(A, B, G, T[j]) * (z * z- I*z);
	testde[5]=-(1.-r*r)*lf[j]*g_intg(A, B, G, T[j]) * (z * z- I*z);
	test[1]=test[0]*r* lf[j]*lf[j]*g_intg(A, B, G, T[j]) * (z * z- I*z);
	test[5]=-test[0]*(1.-r*r)*lf[j]*g_intg(A, B, G, T[j]) * (z * z- I*z);
	for(int i=2; i<5; i++)
	{
		test[i]=0.0;
		testde[i]=0.0;
	}
	for (int k=j; k<tenor_number-1; k++)
	{
		ke[k]=k0;
		sigma[k]=sigma0;
		rho[k]=rho0;
		kappa(r, ke[k], sigma[k], rho[k], lf[j], tenor_number, k, j,Kappa);
		beta(r, sigma[k],lf[j], k, j, b);
		alpha(Kappa,b,rho[k], z,a);
		fd(a,b,z,d);
		fg(a,d,g);
		fA(a,b,d,g, ke[k], sigma[k],j,fa);
		fB(a ,b, d, g, sigma[k], j, fb);
		test[0]=test[0]*exp(fa[0]+fb[0]);
		for(int i=1;i<6; i++)
		{
			testde[i]=testde[i]+fa[i]+fb[i];
		}
	}	
	for (int i=1;i<6; i++)
	{
		test[i]=test[0]*testde[i];
	}	
	
}
void fourier(double r0, double k0, double sigma0, double rho0, double lf0, int tenor_number, int j, double z,complex<double> result[])
{
	complex<double> phi [6];
	complex<double> temp = z*(z-I);
	if((pow(real(temp),2.)+pow(imag(temp),2.))< 1.e-10) 
	{
		result[0]=0.5*(1.-r0*r0)*lf0*lf0*g_intg(A, B, G, T[j]);
		result[1]=-r0*lf0*lf0*g_intg(A, B, G, T[j]);
		result[5]=(1.-r0*r0)*lf0*g_intg(A, B, G, T[j]);
		for (int i=2; i<5; i++)
		{
			result[i] =0.0;
		}
	}
	else
	{
		Phi(r0, k0, sigma0, rho0, lf0, tenor_number, j, z,phi);
		result[0]= (1.0-phi[0])/(z*(z-I));
		for (int i=1; i<6; i++)
		{
			result[i]=-phi[i]/(z*(z-I));
		}
	}
}
void ifft(double r0, double k0, double sigma0, double rho0, double lf0, int tenor_number, int i)
{
	int  Nmax=2*kmax;
	complex<double> ft [6];

	// function input for inverse fourier transform

	// the function value for xi[0]
	xi[0]=0.;
	fourier(r0, k0, sigma0, rho0, lf0, tenor_number, i, xi[0],ft);
	vxi[0]=real(ft[0]);
	vxi1[0]=real(ft[1]);
	vxi2[0]=real(ft[2]);
	vxi3[0]=real(ft[3]);
	vxi4[0]=real(ft[4]);
	vxi5[0]=real(ft[5]);

	// the real part of function value for the last one in xi array
 	xi[kmax]=PI/h;// the last one in xi array
	fourier(r0, k0, sigma0, rho0, lf0, tenor_number, i, xi[kmax],ft);
	vxi[1]=real(ft[0]);
	vxi1[1]=real(ft[1]);
	vxi2[1]=real(ft[2]);
	vxi3[1]=real(ft[3]);
	vxi4[1]=real(ft[4]);
	vxi5[1]=real(ft[5]);
	
	for(int j=1;j<kmax;j++)
	{
		xi[j]=2.*PI*j/h/Nmax;
		fourier(r0, k0, sigma0, rho0, lf0, tenor_number, i, xi[j],ft);

		//real part of input function value for xi[j]
		vxi[2*j] =real(ft[0]);
		vxi1[2*j]=real(ft[1]);
		vxi2[2*j]=real(ft[2]);
		vxi3[2*j]=real(ft[3]);
		vxi4[2*j]=real(ft[4]);
		vxi5[2*j]=real(ft[5]);

		//imaginary part of input function value for xi[j]
		vxi[2*j+1]= imag(ft[0]);
		vxi1[2*j+1]=imag(ft[1]);
		vxi2[2*j+1]=imag(ft[2]);
		vxi3[2*j+1]=imag(ft[3]);
		vxi4[2*j+1]=imag(ft[4]);
		vxi5[2*j+1]=imag(ft[5]);
	}

	realfastfouriertransform(vxi, Nmax, true); // inverse fourier tranform
	realfastfouriertransform(vxi1, Nmax, true); // inverse fourier tranform
	realfastfouriertransform(vxi2, Nmax, true); // inverse fourier tranform
	realfastfouriertransform(vxi3, Nmax, true); // inverse fourier tranform
	realfastfouriertransform(vxi4, Nmax, true); // inverse fourier tranform
	realfastfouriertransform(vxi5, Nmax, true); // inverse fourier tranform

//	the replacement for the output of inverse fourier transform, but not for the direct fourier transform
	for (int n=0; n<kmax; n++)
	{
		ifftresult[kmax-n]=vxi[n]/h;
		ifftresult1[kmax-n]=vxi1[n]/h;
		ifftresult2[kmax-n]=vxi2[n]/h;
		ifftresult3[kmax-n]=vxi3[n]/h;
		ifftresult4[kmax-n]=vxi4[n]/h;
		ifftresult5[kmax-n]=vxi5[n]/h;
	}
	ifftresult[kmax]=0.;
	for (int n=1; n<kmax; n++)
	{
 		ifftresult[n+kmax]=vxi[n]/h;
		ifftresult1[n+kmax]=vxi1[n]/h;
		ifftresult2[n+kmax]=vxi2[n]/h;
		ifftresult3[n+kmax]=vxi3[n]/h;
		ifftresult4[n+kmax]=vxi4[n]/h;
		ifftresult5[n+kmax]=vxi5[n]/h;
	}
}
void calcapprice(double r0, double k0, double sigma0, double rho0, double lf0, int i, int k, double result[] )
{
	ifft( r0, k0,  sigma0, rho0, lf0, TN, i);
	if(strike[k]> forward_rate[i])
	{
		for (int j=kmax; j<2*kmax+1; j++)
		{
			if (log(strike[k]/forward_rate[i])>=(j-kmax)*h && log(strike[k]/forward_rate[i])<=(j+1-kmax)*h)
			{
				result[0]=(T[i+1]-T[i])*zero_bond[i]*forward_rate[i]/(2.*PI)*ifftresult[j+1];
				result[1]=(T[i+1]-T[i])*zero_bond[i]*forward_rate[i]/(2.*PI)*ifftresult1[j+1];
				result[2]=(T[i+1]-T[i])*zero_bond[i]*forward_rate[i]/(2.*PI)*ifftresult2[j+1];
				result[3]=(T[i+1]-T[i])*zero_bond[i]*forward_rate[i]/(2.*PI)*ifftresult3[j+1];
				result[4]=(T[i+1]-T[i])*zero_bond[i]*forward_rate[i]/(2.*PI)*ifftresult4[j+1];
				result[5]=(T[i+1]-T[i])*zero_bond[i]*forward_rate[i]/(2.*PI)*ifftresult5[j+1];
			}
		}
	}
	else 
	{
		for (int j=kmax; j>0; j--)
		{
			if (log(strike[k]/forward_rate[i])<=(j-kmax)*h && log(strike[k]/forward_rate[i])>(j-1-kmax)*h)
			{
				result[0]=(T[i+1]-T[i])*zero_bond[i]*(forward_rate[i]-strike[k])+(T[i+1]-T[i])*zero_bond[i]*forward_rate[i]/(2.*PI)*ifftresult[j-1];
				result[1]=(T[i+1]-T[i])*zero_bond[i]*forward_rate[i]/(2.*PI)*ifftresult1[j-1];
				result[2]=(T[i+1]-T[i])*zero_bond[i]*forward_rate[i]/(2.*PI)*ifftresult2[j-1];
				result[3]=(T[i+1]-T[i])*zero_bond[i]*forward_rate[i]/(2.*PI)*ifftresult3[j-1];
				result[4]=(T[i+1]-T[i])*zero_bond[i]*forward_rate[i]/(2.*PI)*ifftresult4[j-1];
				result[5]=(T[i+1]-T[i])*zero_bond[i]*forward_rate[i]/(2.*PI)*ifftresult5[j-1];
			}
		}
	}
}

void MinimizeObjectiveFunction(double r0, double k0, double sigma0, double rho0, double lf0, int i,double result[])// weighted sum of squares of the corresponding differences between observed maket prices (Black-Scholes prices) and the prices induced by the model, the weights are taken to be proportional to Black-Scholes Vegas: ICV*sqrt(T-t)*exp(-d1*d1/2)/sqrt(2*PI), d1= (ln(ICV/Strike)+sigma*sigma*(T-t)/2)/(sigma*sqrt(T-t))
{
	double capprice [6];
	for(int m=0; m<6; m++)
	{
		capprice[m]=0.0;
		result[m]=0.0;
	}
	for (int k=1; k<SN; k++)
	{
		calcapprice(r0, k0, sigma0,  rho0,  lf0,  i,  k,capprice);
		if (caplet_vola[k][i]*sqrt(T[i])<1.e-10)
		{
			double BlackScholesVega=(log(forward_rate[i]/strike[k])+pow(caplet_vola[k][i],2.0)*T[i]/2.0)/(caplet_vola[k][i]*sqrt(T[i]));
			result[0]=result[0]+(caplet_price1[k][i]-capprice[0])*(caplet_price1[k][i]-capprice[0])/pow(caplet_price1[k][i], 2.0)*BlackScholesVega;
			for (int m=1; m<6; m++)
			{
				result[m]= result[m]-2.*(caplet_price1[k][i]-capprice[0])*capprice[m]/pow(caplet_price1[k][i], 2.0) *BlackScholesVega;
			}
		}
		else 
		{
			double BlackScholesVega=pow(caplet_vola[k][i],2.0)*T[i]/2.0/(caplet_vola[k][i]*sqrt(T[i]));
			result[0]=result[0]+(caplet_price1[k][i]-capprice[0])*(caplet_price1[k][i]-capprice[0])/pow(caplet_price1[k][i], 2.0)*BlackScholesVega;
			for (int m=1; m<6; m++)
			{
				result[m]= result[m]-2.*(caplet_price1[k][i]-capprice[0])*capprice[m] /pow(caplet_price1[k][i], 2.0)*BlackScholesVega;
			}
		}
	}
}

class Functional : public BaseFunctional 
{
	public:
	double evaluateFG(double * x, double * g, int n);
};
double Functional::evaluateFG(double * x, double * g, int n)
{
	double rms;
	double r0 = x[0];
	double k0 = x[1];
 	double sigma0 = x[2];
 	double rho0 = x[3];
	double lf0 = x[4];
	double capprice [6];
	double difference [6];
	MinimizeObjectiveFunction(r0, k0, sigma0, rho0, lf0, maturity, difference);
	rms =difference[0];
	g[0]=difference[1];
	g[1]=difference[2];
	g[2]=difference[3];
	g[3]=difference[4];
	g[4]=difference[5];
	return rms;
}

int main()
{
	TN=42; //maximum number in tenor structure
	SN=7;//number of different strike in the cap-strike martket
	double tenor=0.5; // tenor spacing
	
	//precalibrate result of LMM model, as in (36) of paper
	A=5.001; 
	B=2.000;
	G=2.578;

	allocate_memory(TN, SN);
	initialize(TN, SN);
	InitialCurveAndZeroBond(TN, tenor);
	ReadMarketData(TN, tenor, SN);
	CapletPrice1(TN, SN);
	
	double x[5]= {0.8, 1.40, 1.6900, -0.8000, 3.1};//initial value of parameters
	double l[5]= {0.1, 0.00, 0.0006, -0.9999, 0.1};//lower bound of parameters
	double u[5]= {0.99, 2.00, 2.0000, 0.99999, 4.0};// upper bound of parameters 
	long int ll[5] = {2,2,2,2,2};//number of bounds for each paramter

	for(maturity=TN-2; maturity>0; maturity--)
	{
		if(caplet_price1[1][maturity]!=0.0)
		{
			Functional F;
			double minval =  minimizeBFGS(x,l,u, ll, 5, F, 1.e-2, 1.e-3, 100);
		}
		r0=x[0];
		ke[maturity]=x[1];
		sigma[maturity]=x[2];
		rho[maturity]=x[3];
		lf[maturity]=x[4];
		fprintf(output, "time maturity=%f, r0=%f, ke=%f, sigma=%f, rho=%f, lf=%f\n", maturity*tenor, r0, ke[maturity], sigma[maturity], rho[maturity], lf[maturity]);
		cout<<"time maturity="<<maturity*tenor<<" r="<<r0<<" ke="<<ke[maturity]<<" sigma="<<sigma[maturity]<<" rho="<<rho[maturity]<<" lf="<<lf[maturity]<<endl;
	}
	
	free_memory(TN, SN);
 	fclose (output);
	return 0;
}

