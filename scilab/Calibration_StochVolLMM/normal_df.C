#include "normal_df.h"

 
double gammln(double xx)
{
	double x, y, tmp,ser;
	static double cof[6]={76.1800912947146, -86.50532032941677, 24.01409824083091, -1.2317395722450155, 0.120865097386617e-2, -0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for(j=0; j<5; j++) ser += cof[j]/++y;
	return -tmp+log(2.50662827463100005*ser/x);
}
void gser(double *gamser, double a, double x, double *gln)/* Returns the incomplete gamma function P(a, x) evaluated by its series representation as gamser. Also returns ln\Gamma(a) as gln*/
{
	int n;
	double sum, del, ap;

	*gln=gammln(a);
	if(x<=0.0)
	{
		if(x<0.0) printf("x less that 0 in toutine gser\n");
		*gamser =0.0;
		return;
	}
	else
	{
		ap=a;
		del=sum=1.0/a;
		for(int n=1; n<=ITMAX; n++)
		{
			++ap;
			del *=x/ap;
			sum +=del;
			if(fabs(del) < fabs(sum)*EPS)
			{
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		printf(" a too large, ITMAX too small in routine gser\n");
		return;
	}
}
void gcf(double *gammcf, double a, double x, double *gln) /*Returns the incomplete gamma funcion Q(a, x) evaluated by its continued fraction representation as gammcf. Also returns ln\Gamma(a) as gln.*/
{
	int i;
	double an, b, c, d, del, h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for(i=1; i<ITMAX; i++)
	{
		an= -i*(i-a);
		b+=2.0;
		d=an*d+b;
		if(fabs(d)<FPMIN) d=FPMIN;
		c=b+an/c;
		if(fabs(c)<FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if(fabs(del-1.0)<EPS) break;
	}
	if (i>ITMAX) printf("a too large, ITMAX too small in gcf\n");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
double gammp(double a, double x) /* returns the imcomplete gamma function P(a, x)*/
{
	double gamser, gammcf, gln;
	if (x<0.0|| a<= 0.0) printf("Invalid arguments in routine gammp\n");
	if (x<(a+1.0))
	{
		gser(&gamser, a, x, &gln);
		return gamser;
	}
	else 
	{
		gcf(&gammcf, a, x, &gln);
		return 1.0-gammcf;
	}
}
double gammq(double a, double x) /*Returns the imcomplete function Q(a, x)=1-P(a, x).*/
{
	double gamser, gammcf, gln;
	if (x<0.0|| a<=0.0) printf ("Invalid arguments in routine gammq\n");
	if (x<(a+1.0))
	{
		gser(&gamser, a,x, &gln);
		return 1.0-gamser;
	}
	else
	{
		gcf(&gammcf, a, x, &gln);
		return gammcf;
	}
}
double erff(double x)/* Return the error function erf(x)*/
{
	double gammp(double a, double x);
	return x< 0.0 ? -gammp(0.5, x*x) : gammp(0.5, x*x);
}
double erffc(double x)
{
	double gammp(double a, double x);
	double gammq(double a, double x);
	return x< 0.0 ? 1.0+gammp(0.5, x*x) : gammq(0.5, x*x);
}
double normal_df(double x)
{
	return (0.5+0.5*erff(x/sqrt(2)));
}


