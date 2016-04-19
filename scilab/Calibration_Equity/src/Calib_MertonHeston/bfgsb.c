#include "bfgsb.h"

extern  void setulb_(int*, int*, double*, double*, double*, int*, double*, double*, double*, double*,
					double*, int*, char*, int*, char*, int*, int*, double*);
extern  void inittask_(char*);
extern  void inittask2_(char*);

static int minix = 0;


int bfgsb(int dimx ,double *x,int *nbd,double *xmin,double *xmax)
{
  int i,dimx1,maxiter;
  double fx,grad[dimx];
  int mg,nbwa,niwa,*iwa,iprint;
  double factr,pgtol,*wa;
  char task[60], csave[60];
  int lsave[4],isave[44];
  double dsave[29];
  double *bestx,bestfx;
  

  bestx = (double*)malloc(dimx*sizeof(double));

  if(minix == 1)
	{
	  fx = costFunction1(dimx,x);
	  gradFunction1(dimx,x,grad);
	} else {
	  fx = costFunction(dimx,x);
	  gradFunction(dimx,x,grad);
	}
  
  //for(i=0;i<dimx;i++) printf("nbd[%d]=%d \n",i,nbd[i]);
  
  for(i=0;i<dimx;i++) bestx[i] = x[i];
  bestfx = fx;
  
  //=================
  // Init des parametres du BFGS a memoire limitee
  //
  dimx1 = dimx;
  // number of corrections used in the limited memory matrix
  mg = mg0; 
  // ~ tolerance sur la variation relative (* presicion machine )de f
  //: e12 faaible precision, e7 moyenne, et e1 tres grande
  factr = factr0; 
  // tolerace sur la norme du gradient
  pgtol = pgtol0;
  // tableaux temporaires
  nbwa = (2*mg + 4)*dimx + 12*mg*mg + 12*mg ;
  niwa = 3*dimx ;
  wa = (double *) malloc(nbwa*sizeof(double));  
  iwa = (int *) malloc(niwa*sizeof(int));  
  // pour le type d'affichage a chaque iteration
  iprint = iprint0;
  // nombre maximum d'iterations
  maxiter = maxiter0;
  //=================
  inittask_(task); 
  setulb_(&dimx1, &mg, x, xmin, xmax, nbd, &fx, grad, &factr, &pgtol,
		  wa, iwa, task, &iprint, csave, lsave, isave, dsave);



  
  while (memcmp(task,"FG",2) == 0 || memcmp(task,"NEW_X",5) == 0)
	{
	  if (memcmp(task,"FG",2) == 0) // si les 2 premiers caracteres de task sont FG
	    {
		  //
		  if(minix == 1)
			{
			  fx = costFunction1(dimx,x);
			  gradFunction1(dimx,x,grad);
			} else {
			  fx = costFunction(dimx,x);
			  gradFunction(dimx,x,grad);
			}
		  // ATTENTION la ligne qui suit est une vrai Connerie
		  //grad[1] = 10.*grad[1];
		  //
	    }
	  if (memcmp(task,"NEW_X",5) == 0 && isave[29]>=maxiter)
	    //isave[29] en C ou isave(30) en Fortran correspond au nb d'iter. courant
	    {
	      inittask2_(task); // task = 'STOP: Maximum number of iterations reached'
	    }
	  setulb_(&dimx1, &mg, x, xmin, xmax, nbd, &fx, grad, &factr, &pgtol,
			  wa, iwa, task, &iprint, csave, lsave, isave, dsave);
	  if(fx<bestfx)
		{
		  for(i=0;i<dimx;i++) bestx[i] = x[i];
		  bestfx = fx;
		}

	  //printf("Dans bfgsb.c : fx = %f \n",fx);
	  
	}
  

  for(i=0;i<dimx;i++) x[i] = bestx[i];
  //
  free(wa);
  free(iwa);
  free(bestx);
  //
  return 0;
  
}

void initbfgsb(int _minix)
{
  minix = _minix;  
}

