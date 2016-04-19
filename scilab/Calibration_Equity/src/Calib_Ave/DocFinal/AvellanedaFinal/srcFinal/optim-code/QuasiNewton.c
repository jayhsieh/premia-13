//
// MATHFI Project, Inria Rocquencourt.
// Vincent Barette, June 2002.
//


/* 
** FILE
** *****************************************************************
** 
** QuasiNewton.c
** 
** 
** PURPOSE
** 
** Implements a Quasi-Newton unconstrained optimizer. Uses a BFGS update
** of the Cholesky factors of the approximation of the hessian. An initial
** scaling is made in the first iteration.
** 
** 
*****************************************************************
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "QuasiNewton.h"
#include "lineSearch.h"
#include "BFGSupdate.h"
#include "cholesky.h"
#include "stopping.h"
#include "testing.h"



FILE *fich;



/* 
** FUNCTION
** *****************************************************************
** 
** macheps()
** 
** 
** DESCRIPTION
** 
** Computes the machine epsilon. Return a double used by the update
** of cholesky factors.
** 
** 
** 
*****************************************************************
*/


double macheps()
{
  double e;
  e=1.;
  while ( (1+e)!=1 )
    {
      e=e/2.;
    }
  return e;
}


/* 
** FUNCTIONS
** *****************************************************************
** 
** printPointGradC, printVect, printFinalPoint, saveInFile
** 
** 
** DESCRIPTION
** 
** Several functions printing some information in the terminal during the
** algorithm, or saving values in the file "data.out".
** 
** 
** 
*****************************************************************
*/

void printPointGradC(int n, int counter, int verbosity, double* x0,double* g0,struct chol *C)
{
  int i;
  if (verbosity>=1 ) {printf("\nITERATION %d:",counter);}
  if (n<=4 && verbosity==3)
    {for (i=0;i<n;i++){printf(" %e  ",x0[i]);}}
  if (verbosity>=1 ) {printf("\n");}
  if (n<=4 && verbosity==3 )
    {printf("Gradient g%d= ",counter);
     for (i=0;i<n;i++){printf(" %e  ",g0[i]);}
     printf("\n");}
  if ( n<=3 && verbosity==3 ){printCholM(C);}
  if (n<=4 && verbosity==3 )
    {printf("Elements of the D matrix= ");
     for (i=1;i<=n;i++) {printf(" %e  ",elementDi(i,C));}
     printf("\n");}
}

void printVect(char* c, int n, int verbosity, int counter, double* vect)
{
  int i;
  if (n<=4 && verbosity==3)
    {printf(c,counter);
     for (i=0;i<n;i++) {printf(" %e  ",vect[i]);}
     printf("\n");}
}

void printFinalPoint(int n, int verbosity, int counter, double* x0)
{ 
  int i;
  if (n<=4 && verbosity==3)
    {printf("\nFINAL POINT x%d:",counter);
     for (i=0;i<n;i++){printf(" %e  ",x0[i]);}}
  printf("\n");
}

void printStoppingCriteria(double relativeGradient, double relativeChangeIn_x, int verbosity, double gradtol, double steptol)
{
  if (verbosity>=3)
    {printf("Stopping criteria:\n");
     printf("\trelative gradient    (<%e ?): %e\n", gradtol, relativeGradient);
     printf("\trelative change in x (<%e ?): %e\n", steptol, relativeChangeIn_x);}
}

void printNormGradient(int n, double* g, double N2g0, int verbosity)
{
  double N2g=0;
  double d=0;
  double f;
  int i;

  for (i=0;i<n;i++) { N2g+= g[i]*g[i];}
  N2g=sqrt(N2g);
  f=N2g/N2g0;
  if (verbosity>=2) printf("||grad||/||grad0||= %e\n",f);
}

void saveInFile(int n, FILE *fich, double* x0)
{
  int i;
  for (i=0;i<n;i++) {fprintf(fich,"%lf\n",x0[i]);}
}


/* 
** FUNCTION
** *****************************************************************
** 
** scaling
** 
** 
** DESCRIPTION
** 
** Used in the first iteration of the Quasi-Newton algorithm. Multiplies
** the hessian matrix by a scaling factor.
** 
** 
** 
*****************************************************************
*/

void scaling(int n,double* x0,double* x1,double* g0,double* g1,struct chol *C,int verbosity)
{
  int i;
  double ys=0;
  double ss=0;
  double mu;
  double si;
  for(i=0;i<n;i++)
    {
      si = x1[i]-x0[i];
      ys+= (g1[i]-g0[i])*si;
      ss+= si*si;
    }
  mu=ys/ss;
  if (verbosity>=3) {printf("Scaling factor= %e\n",mu);}
  scalarMult(C,mu);
}


/* 
** FUNCTION
** *****************************************************************
** 
** void QuasiNewton(int n, 
**		    double* starting_x, 
**		    double (*costFunction)(double*), 
**		    void  (*gradCostFunction)(double*,double*),
**		    double* sol,
**		    double gradtol,
**		    double steptol,
**		    int maxCounter,
**		    int verbosity
**		    )
** 
** 
** DESCRIPTION
** 
** Minimize the cost-function given the starting point "starting_x".
** The optimal point is stocked in "sol".
** 
** The method used is that of Quasi-Newton BFGS update.
** It is terminated by either reaching the maximum number of iterations or
** satisfying one of the two tolerance criteria:
**    - gradtol : on the relative gradient
**    - steptol : on the relative change in successive values of x.
**
** Inputs :
**      int n: size of the problem.
**	starting_x: initial point.
**	costFunction : the cost-funtion to minimize.
**	gradCostFunction(x,g): stocks in g the gradient of the cost-function at x.
**	sol: receives the final solution.
**	gradtol: tolerance on the relative gradient.
**	steptol: tolerance on the relative change of x.
**	maxCounter: maximum number of iterations.
**	verbosity: level of information choosen by the user.
**      saveSuccessiveXinFile: save succesive x0 in the file "data.out".
**
** 
*****************************************************************
*/

void QuasiNewton(int n,
		 double* starting_x, 
		 double (*costFunction)(double*), 
		 void  (*gradCostFunction)(double*,double*),
		 double* sol,
		 double gradtol,
		 double steptol,
		 int maxCounter,
		 int verbosity,
		 int saveSuccessiveXinFile
		 )
{
  
  int counter=0;
  int totalCallsJ=0;
  int totalCallsGradJ=0;
  int endQuasiNewton=0;
  int endGradient=0;
  int endStepX=0;
  int i,j;

  double epsilon=macheps();
  double typf=0;
  double relativeGradient;
  double relativeChangeIn_x;
  double alpha;//the step given by the line search
  double fx0;
  double finit;
  double N2g0=0;

  struct chol *C;

  double* x0;
  double* x1;
  double* g0;
  double* g1;
  double* v;
  double* searchdir;//the direction of the line search
  double* typx;
  double scalaire=0.0,N2searchdir=0.0,N2g1=0.0;

	if (saveSuccessiveXinFile==1) {fich=fopen("data.out","w");}

  C = (struct chol *) malloc(n*sizeof(struct chol));
  C->size = n;
  C->l = (double *) malloc( (n*(n-1))/2 * sizeof(double));
  C->d = (double *) malloc( n * sizeof(double));
  init(C);

  x0        = (double *) malloc(n*sizeof(double));
  x1        = (double *) malloc(n*sizeof(double));
  g0        = (double *) malloc(n*sizeof(double));
  g1        = (double *) malloc(n*sizeof(double));
  v         = (double *) malloc(n*sizeof(double));
  searchdir = (double *) malloc(n*sizeof(double));
  typx      = (double *) malloc(n*sizeof(double));

  for (j=0;j<n;j++)
    {
      x0[j]=0;
      x1[j]=0;
      g0[j]=0;
      g1[j]=0;
      v[j]=0;
      searchdir[j]=0;
    }

  for (j=0;j<n;j++) { x0[j]=starting_x[j]; }

  gradCostFunction(x0,g0);

  // computes N2(g0)
  for (j=0;j<n;j++) { N2g0 += g0[j]*g0[j];}
  N2g0=sqrt(N2g0);

  //typical values of x[i] components. Used in the relative gradient.
  for (j=0;j<n;j++) { if (x0[j]!=0) {typx[j]=x0[j];} else {typx[j]=1;} }


  /*********************************************************************************/
  /********************** beginning of the while loop  *****************************/

  while ( counter < maxCounter+1 && endQuasiNewton==0 )
    {

      //prints some information
      printPointGradC(n,counter,verbosity,x0,g0,C);

      //save successive x0 in a file
      if (saveSuccessiveXinFile==1) {saveInFile(n,fich,x0);}

      //computes the search direction solving "H*d=-grad"
      cholSolveOpp(C,searchdir,g0);

      //prints some information
      printVect("Direction d%d= ",n,verbosity,counter,searchdir);

      //the line-search. Return the step alpha.
      alpha=Wolfe_cubic( n,costFunction,gradCostFunction,x0,searchdir,&fx0,&totalCallsJ,&totalCallsGradJ,verbosity);

      if (counter==0) {finit=fx0;}
      //printf("alpha %f\n",alpha);
      //The next point x1 = x0 + step*direction
      for (i=0;i<n;i++) { x1[i] = x0[i] + alpha*searchdir[i];}

      //prints some information
      printVect("point x%d= ",n,verbosity,counter+1,x1);

      //A typical value "typf" of the cost-function. Used in the relative gradient.
      if (counter==0) { if (fx0!=0) {typf=fx0;} else {typf=1;} }

      //Relative Gradient stopping criteria
      relativeGradient=relgrad(n,x0,g0,fx0,typx,typf);

      //Relative change in x  stopping criteria.
      relativeChangeIn_x = relx(n,x1,x0,typx);

      //prints some information
      printStoppingCriteria(relativeGradient, relativeChangeIn_x, verbosity,gradtol,steptol); printNormGradient(n, g0, N2g0,verbosity);

	/*		//	AJOUT DE marouen : affichage de quelques stat
			scalaire=0.0;
			N2g1=0.0;
			N2searchdir=0.0;
			for (j=0;j<n;j++) {  scalaire+= g0[j]*searchdir[j];}
			for (j=0;j<n;j++) {  N2g1+= g0[j]*g0[j];N2searchdir+= searchdir[j]*searchdir[j];}
			N2searchdir=sqrt(N2searchdir);
			N2g1=sqrt(N2g1);
			scalaire =scalaire/(N2searchdir*N2g1);
      //ajout X pour calcul des erreurs relatives a chaque iterations

			erreurPrixOptim=0.0;
      for( k = 0 ; k<SOption.size ;k++)
			{
				prixOptim=optimPrice(SOption.strik[k],param.S0,SOption.maturite[k],param, flagOption);
				//tempVol = ImpliedVol(param.S0, param.r, param.dividende, flag , 0 , SOption.strik[k], SOption.maturite[k], prixOptim,0.00000001);
				erreur =fabs(prixOptim-SOption.prix[k])/SOption.prix[k];
				//printf("prix initial %f prix %f erreur %f impVol %f \n",SOption.prix[i] ,prixOptim, erreur, tempVol );
				//printf("%f %f %f %f %f %f\n",SOption.strik[k],SOption.maturite[k],SOption.prix[k],prixOptim,erreur,tempVol);
				erreurPrixOptim +=erreur;
			}

			//fin ajout X */




      //Test of stopping criteria
      if ( relativeGradient < gradtol )
				{endGradient=1; counter++; break;}
      if ( relativeChangeIn_x < steptol) 
				{for (i=0;i<n;i++) { x0[i] = x1[i];}
	 		endStepX=1;counter++;break;}
			//if (counter==6){endGradient=1; counter++; break;}
      //g1 receives the gradient of the cost-function at x1
      gradCostFunction(x1,g1);
			
      //scaling for the first iteration.
      if (counter==0) { scaling(n,x0,x1,g0,g1,C,verbosity); }

      //BFGS update
      if (alpha!=0) {BFGSupdate(C,x1,x0,g1,g0,n,epsilon,&endQuasiNewton,verbosity);}

      //x0 receives x1 and g0 receives g1
      for (i=0;i<n;i++)
			{
			  x0[i] = x1[i];
	  		g0[i] = g1[i];
			}
      
      counter++;

			
    }
  
  /************************ end  of the while loop  *******************************/  
  /*********************************************************************************/
  

  //sol receives the optimal point
  for (i=0;i<n;i++) { sol[i] = x0[i]; }
  
  //prints some final information
  if (verbosity>=1 ) { 
    printf("\n***********************END*********************************\n");
    printf("#iterations realized : %d",counter-1);
    printf("\nStopping criteria    : ");
    if (counter==maxCounter+1) {printf("Maximum number of iterations reached");}
    if (endGradient==1) {printf("relative gradient < %e",gradtol);}
    if (endStepX==1) {printf("relative change in x < %e",steptol);}
    printf("\nDecrease of J        : ");
    printf("%e --> %e",finit,fx0);
    printf("\n#calls               : J --> %d  gradJ --> %d\n",totalCallsJ,totalCallsGradJ);
    printf("***********************************************************\n");

    printFinalPoint(n, verbosity, counter, x0);
  }
  
  //memory disallocation
  free(searchdir);
  free(x0);
  free(x1);
  free(g0);
  free(g1);
  free(v);
  free(typx);
  free(C->l);
  free(C->d); 
  free(C);
	
  if (saveSuccessiveXinFile==1) {fclose(fich);}
}



/* 
** FUNCTION
** *****************************************************************
** 
** testQuasiNewton(testProblem T, int verbosity)
** 
** 
** DESCRIPTION
** 
** tests teh optimizer with classical problems (the structure testProblem
** is defined in testing.h)
** 
** 
** 
*****************************************************************
*/


void testQuasiNewton(testProblem T, int verbosity)
{
  double* sol = (double *) malloc(T.size*sizeof(double));
  QuasiNewton(T.size, T.startingPoint, T.costFunction, T.gradCostFunction, sol, 1e-6,1e-6,200,verbosity,1);
  free(sol);
}





