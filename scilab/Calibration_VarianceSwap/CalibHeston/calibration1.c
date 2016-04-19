
#include "utils/optim.h"

extern double costFunction(int dimx, double *x);
extern void  gradFunction(int dimx,double *x, double *grad);
extern double costFunctionVS(int dimx, double *x);
extern void  gradFunctionVS(int dimx,double *x, double *grad); 
extern int TypeModel;
extern int ifChangeVar;

int checkRanges(int dimx, double *x)
{
	int ifOK, i;

	ifOK = 1;
	for(i=0; i<4; i++) 
	{
		if(x[i]<0.0)
		{
			x[i]=1e-04;
			ifOK = 0;
		}
	} 
	if(fabs(x[4])>0.9999)
		{
			x[4] = fabs(x[4])/x[4]*0.9;
			ifOK = 0;
		}
	return ifOK;
}
//=========================== OPTIM CALIBRATION ============================
void optim_calibration(int *ind, int *n, double *x, double *f, double *g, opt_simul_data *optim_data)
{
  double grad[5];
  double xorig[5];
  int i;

  if( (TypeModel == 111)||(!ifChangeVar)) 
  {
	if(!checkRanges(*n, x)) {printf("Some of parameters are out of range!\n");} 
	for(i=0; i<5; i++) xorig[i]=x[i];
  }
  else if(ifChangeVar)
  { // if  variables are changed 
	xorig[0] = 0.5 + atan(x[0])/M_PI;
	xorig[1]  = 30.0*(0.5 + atan(x[1])/M_PI);
	xorig[2]  = 0.5 + atan(x[2])/M_PI;
	xorig[3] = 0.5 + atan(x[3])/M_PI;
	xorig[4]    = 2.0*atan(x[4])/M_PI;
  }

  printf("current params %d %f %f %f %f %f\n", TypeModel, xorig[0], xorig[1], xorig[2], xorig[3], xorig[4]);
  
  if (*ind == 2 || *ind == 4 )
    {
      /* cost */
       if(TypeModel==111) { *f  = costFunctionVS(*n,x);}
	else {	
		*f  = costFunction(*n, xorig);
             }
    }
   
  if (*ind == 3 || *ind == 4 )
    {
      /* gradient */
	if(TypeModel==111) { gradFunctionVS(*n,x,grad); 
	g[0] = grad[0];
      g[1] = grad[1];
      g[2] = grad[2];
      g[3] = grad[3];
      g[4] = grad[4];  }
	else {
		gradFunction(*n,xorig,grad);
		g[0] = grad[0];
		g[1] = grad[1];
		g[2] = grad[2];
		g[3] = grad[3];
		g[4] = grad[4];
		if(ifChangeVar)
		{
			g[0] = g[0]/M_PI/(1.0+x[0]*x[0]);
			g[1] = g[1]*30.0/M_PI/(1.0+x[1]*x[1]);
			g[2] = g[2]/M_PI/(1.0+x[2]*x[2]);;
			g[3] = g[3]/M_PI/(1.0+x[3]*x[3]);
			g[4] = g[4]*2.0/M_PI/(1.0+x[4]*x[4]);
		}
		if(TypeModel==11) {// to fix parameters fitted by Variance Swaps
			g[0] =0.0;
			g[1] =0.0;
			g[2] =0.0;
			}
		}
	}
}

//=========================== OPTIM SCILAB =========================  
int optim_scilab( int n, double x[])
{
  int mem=10, nwork, niter=100, nsim=500, impress=20, io=0, mode=0, quatre=4;
  
  double f,g[5], dxmin=1.e-9,df1=1.0,epsg=1.e-12;
   double *work;  
  
  /* working arrays  */
  nwork = (4*n + mem*(2*n + 1));
  if ((work = malloc(nwork*sizeof(double))) == NULL) return 1;
 
 optim_calibration(&quatre,&n,x,&f,g,NULL);

  optim_n1qn3 (optim_calibration,optim_fuclid,optim_ctonb,optim_ctcab,&n,x,&f,g,&dxmin,&df1,&epsg,&impress,&io,&mode,&niter,&nsim,work,&nwork,NULL);
  
  return 0;
}

