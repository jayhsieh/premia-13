
#include "spline.h"

double cubspline(int Metode, double xi, double *xx,double *yy, int N)
{
	int i;
	double yi;
	double *x;
	double *y;
	double *y2;
	int j;
	
	if( Metode == 1 ) 
	{
		//Numerical Recipes are 1 based
		j = 0;
	}
	else
	{	//Others are 0 based
		j = 0;
	}
	
	x=(double *) malloc(N*sizeof(double));
	y=(double *) malloc(N*sizeof(double));
	for( i = 0 ; i<N ;i++)
	{
		//if( strcmp((char*)&yy[i] ,"")==0)
		//if( yy[i]!=0)
		{
			j = j + 1;
			x[j-1] = xx[i];
			y[j-1] = yy[i];
		}
	}
	
	if( Metode == 1)
	{
		//NR cubic spline
		//Get y2
		y2=(double *) malloc(N*sizeof(double));
		spline(x, y, N, pow(10,30), pow(10,30), y2);
		//Get y
		splint(x, y, y2, N, xi, &yi);
	}
	else if( Metode == 3)
	{
		//Own cubic spline
		yi = SplineX3(xi, x, y,N);
	}
	free(x);
	x=NULL;
	free(y);
	y=NULL;
	if( Metode == 1)
	{
		free(y2);
		y2=NULL;
	}
	
	return yi;
}



void spline(double *x, double *y, int N,double  yp1, double ypn ,double *y2)
{
	
	//Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., y i = f(xi), with
	//x1<x2< :::<xN , and given values yp1 and ypn for the first derivative of the inter-
	//polating function at points 1 and n, respectively, this routine returns an array y2(1:n) of
	//length n which contains the second derivatives of the interpolating function at the tabulated
	//points xi. If yp1 and/or ypn are equal to 1 * 10^30 or larger, the routine is signaled to set
	//the corresponding boundary condition for a natural spline, with zero second derivative on
	//that boundary.
	//Parameter: NMAX is the largest anticipated value of n.
	
	int Nmax;
	int i;
	int k ;
	double p;
	double qn;
	double sig;
	double un;
	double  *u;
	
	Nmax = 500;
	u=(double*)malloc(Nmax*sizeof(double));
	//The lower boundary condition is set either to be natural
	if (yp1 > 9.9E+29)
	{
        y2[0] = 0;
        u[0] = 0;
	}
	else
	{
        //or else to have a specicied first derivative.
        y2[0] = -0.5;
        u[0] = (3 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
	}
	
	//This is the decomposition loop of the tridiagonal
	//algorithm. y2 and u are used for temporary
	//storage of the decomposed factors.
	
	for (i = 1; i< N - 1;i++)
	{
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
		p = sig * y2[i - 1] + 2;
        y2[i] = (sig - 1) / p;
        u[i] = (6 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1])) / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
	}
	
	//The upper boundary condition is set either to be natural
	if (ypn > 9.9E+29)
	{
        qn = 0;
        un = 0;
	}
	else
	{
        //or else to have a specified first derivative.
        qn = 0.5;
        un = (3 / (x[N-1] - x[N - 2])) * (ypn - (y[N-1] - y[N - 2]) / (x[N-1] - x[N - 2]));
	}
	y2[N-1] = (un - qn * u[N - 2]) / (qn * y2[N - 2] + 1);
	
	//This is the backsubstitution loop of the tridiagonal algorithm.
	for( k = N - 2;k>=0 ;k--)
	{
        y2[k] = y2[k] * y2[k + 1] + u[k];
	}
	free(u);
}
      
void splint(double *xa, double *ya,double  *y2a,int N ,double x, double *y )
{
	
	//Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
	//xai //s in order), and given the array y2a(1:n), which is the output from spline above,
	//and given a value of x, this routine returns a cubic-spline interpolated value y.
	
	int k;
	int khi;
	int klo;
	
	double A;
	double B;
	double h;
	
	//We will the right place in the table by means of bisection.
	klo = 1;
	khi = N;
	
	while(khi - klo > 1)
	{
        k = (khi + klo) / 2;
		if (xa[k-1] > x)
		{
			khi = k;
		}
        else
		{
			klo = k;
		}
	}
	
	//klo and khi now bracket the input value of x.
	h = xa[khi-1] - xa[klo-1];
	if (h == 0) printf("bad xa input in splint\n");
	
	//Cubic spline polynomial is now evaluated.
	A = (xa[khi-1] - x) / h;
	B = (x - xa[klo-1]) / h;
	*y = A * ya[klo-1] + B * ya[khi-1] + ((pow(A,3) - A) * y2a[klo-1] + (pow(B,3) - B) * y2a[khi-1]) * (pow(h , 2) ) / 6;
	
}

double SplineX3(double x , double *xx, double *yy, int N)
{
	//|-------------------------------------------------------------------------------
	//| Function returns y value for a corresponding x value, based on cubic spline.
	//| Will never oscillates or overshoot. No need to solve matrix.
	//| Also calculate constants for cubic in case needed (for integration).
	//|
	//| xx(0 to No_of_lines) is x values
	//|    * Must be unique (no two consequetive ones the same)
	//|    * Must be in ascending order
	//|    * No of lines = Number of points - 1
	//| yy(0 to No_of_lines) is y values
	//|
	//| Uses function dxx to prevent div by zero.
	//|
	//| Developer: C Kruger, Guildford, UK
	//| Date: December 2001
	//|-------------------------------------------------------------------------------
	
	int i;
	int j;
	int Nmax;
	int Num;
	
	//1st and 2nd derivative for left and right ends of line
	double *gxx;
	double *ggxx;
	
	//Constants for cubic equations
	double A ;     //Also for linear extrapolation
	double B ;     //Also for linear extrapolation
	double C ;
	double D ;
	
	//Number of lines = points - 1
	Nmax = N-1;
	gxx=(double *)malloc(2*sizeof(double));
	ggxx=(double *)malloc(2*sizeof(double));
	
	//(1a) Find LineNumber or segment. Linear extrapolate if outside range.
	Num = 0;
	if( x < xx[0] || x > xx[Nmax])
	{
		//X outisde range. Linear interpolate
		//Below min or max?
		if (x < xx[0])  Num = 1;
		else Num = Nmax;
		B = (yy[Num] - yy[Num - 1]) / dxx(xx[Num], xx[Num - 1]);
		A = yy[Num] - B * xx[Num];
		return  A + B * x;
	}
	
	//(1b) Find LineNumber or segment. Linear extrapolate if outside range.
	else
	{
		//X in range. Get line.
		for( i = 1;i<Nmax+1;i++)
		{
			if( x <= xx[i])
			{
				Num = i;
				break;
			}
		}
	}
	
	//(2) Calc first derivative (slope) for intermediate points
	for( j = 0;j<= 1;j++)          //Two points around line
	{
		i = Num - 1 + j;
		if( i == 0 || i == Nmax)
		{
			//Set very large slope at ends
			gxx[j] = pow(10 , 30);
		}
		else if ((yy[i + 1] - yy[i] == 0) || (yy[i] - yy[i - 1]== 0) )
		{
			//Only check for 0 dy. dx assumed NEVER equals 0 !
			gxx[j] = 0;
		}
		else if( ((xx[i + 1] - xx[i]) / (yy[i + 1] - yy[i]) + (xx[i] - xx[i - 1]) / (yy[i] - yy[i - 1])) == 0) 
		{
			//Pos PLUS neg slope is 0. Prevent div by zero.
			gxx[j]= 0;
		}
		else if( (yy[i + 1] - yy[i]) * (yy[i] - yy[i - 1]) < 0 )
		{
			//Pos AND neg slope, assume slope = 0 to prevent overshoot
			gxx[j]= 0;
		}
		else
		{ //Calculate an average slope for point based on connecting lines
			gxx[j] = 2.0 / (dxx(xx[i + 1], xx[i]) / (yy[i + 1] - yy[i]) + dxx(xx[i], xx[i - 1]) / (yy[i] - yy[i - 1]));
		}
	}
	
	//(3) Reset first derivative (slope) at first and last point
	if (Num == 1)
	{
		//First point has 0 2nd derivative
		gxx[0] = 3.0 / 2 * (yy[Num] - yy[Num - 1]) / dxx(xx[Num], xx[Num - 1]) - gxx[1] / 2.0;
	}
	if( Num == Nmax)
	{
		//Last point has 0 2nd derivative
		gxx[1] = 3.0 / 2 * (yy[Num] - yy[Num - 1]) / dxx(xx[Num], xx[Num - 1]) - gxx[0] / 2.0;
	}
	
	//(4) Calc second derivative at points
	ggxx[0] = -2.0 * (gxx[1] + 2 * gxx[0]) / dxx(xx[Num], xx[Num - 1]) + 6 * (yy[Num] - yy[Num - 1]) / pow(dxx(xx[Num], xx[Num - 1]),2);
	ggxx[1] = 2.0 * (2 * gxx[1] + gxx[0]) / dxx(xx[Num], xx[Num - 1]) - 6 * (yy[Num] - yy[Num - 1]) / pow(dxx(xx[Num], xx[Num - 1]),2);
	
	//(5) Calc constants for cubic
	D = 1.0 / 6 * (ggxx[1] - ggxx[0]) / dxx(xx[Num], xx[Num - 1]);
	C = 1.0 / 2 * (xx[Num] * ggxx[0] - xx[Num - 1] * ggxx[1]) / dxx(xx[Num], xx[Num - 1]);
	B = (yy[Num] - yy[Num - 1] - C * (pow(xx[Num], 2) - pow(xx[Num - 1], 2)) - D * (pow(xx[Num],3) - pow(xx[Num - 1],3))) / dxx(xx[Num], xx[Num - 1]);
	A = yy[Num - 1] - B * xx[Num - 1] - C * pow(xx[Num - 1],2) - D * pow(xx[Num - 1],3);
	
	free(gxx);
	free(ggxx);
	//Return function
	return A + B * x + C * pow(x,2) + D * pow(x,3);
	
	////Alternative method based on Numerical Recipes.
	////Shorter but does not calc cubic constants A, B, C, D
	//i = Num
	//A = (xx(i) - x) / (xx(i) - xx(i - 1))
	//B = 1 - A
	//Cy = 1 / 6 * (A ^ 3 - A) * (6 * (yy(i) - yy(i - 1)) - 2 * (gxx(i) + 2 * gxx(i - 1)) * (xx(i) - xx(i - 1)))
	//Dy = 1 / 6 * (B ^ 3 - B) * (2 * (2 * gxx(i) + gxx(i - 1)) * (xx(i) - xx(i - 1)) - 6 * (yy(i) - yy(i - 1)))
	////Return function
	//SplineX3 = A * yy(i - 1) + B * yy(i) + Cy + Dy
	
}

double dxx(double x1, double x0)
{
	//Calc Xi - Xi-1 to prevent div by zero
	
	double dxx = x1 - x0;
	if (dxx == 0)  dxx = pow(10 , 30);
	return dxx;
}

int Sort(double *x,double*y, int size)
{
    int i, j, inc;
	
    double *tabx,*taby;
	double * array;
	double *w;
    int n = size;
    int nCols = 2;
	double v;
    
	
	tabx =(double*) malloc ( size*sizeof(double));
	for(i=0;i<size;i++)
		tabx[i] = x[i];
	
	taby =(double*) malloc ( size*sizeof(double));
	for(i=0;i<size;i++)
		taby[i] = y[i];
    
	array =(double*) malloc ( size*sizeof(double));
	for(i=0;i<size;i++)
		array[i] = x[i];
	
	//A_Vector w(nCols, 0.);
	w =(double*) malloc ( nCols*sizeof(double));
	for(i=0;i<nCols;i++)
		w[i] = 0.0;
	
	
	
    inc = 1;
	
    do
    {
        inc *= 3;
        inc++;
    }
    while (inc <= n);
	
    do
    {
        inc /= 3;
		
        for (i=inc+1; i<=n; i++)
        {
            v = array[i-1];
			
            w[0] = tabx[i-1];
			w[1] = taby[i-1];
            
			j=i;
			
            while (array[j-inc-1] > v)
            {
                array[j-1] = array[j-inc-1];
				
				tabx[j-1] = tabx[j-inc-1];
				taby[j-1] = taby[j-inc-1];
                
                j -= inc;
				
                if (j <= inc) break;
            }
			
            array[j-1] = v;
			tabx[j-1] = w[0];
			taby[j-1] = w[1];
			
        }
    }
    while (inc>1);
	
	for(i=0;i<size;i++)
	{
		x[i] = tabx[i];
		y[i] = taby[i];
	}
    
	free(w);
	free(tabx);
	free(taby);
	free(array);
	
    return(1);
}
