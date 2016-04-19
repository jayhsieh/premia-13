
#include <fstream>
#include <iomanip>
#include <cstring>
#include <iostream>
#include "min.h"
#include "f2c.h"


using namespace std;

double GradTester::evaluateFG(double * x, double * g, int n){
	double * g2 = new double[n];
	double r;
	for(int i=0; i<n; i++){
		cout << '.';
		x[i] += step;
		r = bf.evaluateFG(x,g2,n);
		x[i]-=2.0*step;
		g[i] = (r-bf.evaluateFG(x,g2,n))/2.0/step;
		x[i]+=step;
	}
	cout << endl;
	r = bf.evaluateFG(x,g2,n);
	double sum=0;
	for(int i=0; i<n; i++) sum+=(g2[i]-g[i])*(g2[i]-g[i]);
	if(true){
		cerr << "GradTester warning: difference of gradients too large."<<endl;
		cerr << "Gradients dumped to file gradtester.txt"<<endl;
		ofstream ofile("gradtester.txt");
		ofile << "%" << setw(10) << "X"  << setw(10) <<"g"<<setw(10) <<"gnum"<<endl;
		for(int i=0; i<n; i++) ofile << setw(2)<< "" << setw(9)<<x[i] << " "<<setw(9)<< g2[i]<< " "<<setw(9)<< g[i]<<endl;
		throw 1;
	}
	delete[] g2;
	return r;
}


/*extern "C" int setulb_(integer *n, integer *m, doublereal *x, 
	doublereal *l, doublereal *u, integer *nbd, doublereal *f, doublereal 
	*g, doublereal *factr, doublereal *pgtol, doublereal *wa, integer *
	iwa, char *task, integer *iprint, char *csave, logical *lsave, 
	integer *isave, doublereal *dsave, ftnlen task_len, ftnlen csave_len);
*/
extern "C" int setulb_(long int *n, long int *m, double *x, 
	double *l, double *u, long int *nbd, double *f, double 
	*g, double *factr, double *pgtol, double *wa, long int *
	iwa, char *task, long int *iprint, char *csave, long int  *lsave, 
	long int *isave, double *dsave, short task_len, short csave_len);

double minimizeBFGS(double * x, double * l, double * u, long int * nbd, int n, BaseFunctional & func, double factr, double pgtol, int maxit){

// Limited Memory BFGS with linear constraints
// Shell function for the FORTRAN routine
// x contains the initial guess
// l contains lower bounds; u contains upper bounds
// nbd describes the nature of bounds
// 0 ==== no bounds
// 1 ==== lower bound
// 2 ==== lower and upper 
// 3 ==== only upper

    char task[60];
    char csave[60];
    logical lsave[4];
    const int mmax =17;
    integer  iprint;
	integer m;
    integer * iwa = new integer[3*n];
    integer isave[44];
    double f;
    double  * g = new double[n];
    double dsave[29];
    double * wa = new double [2*mmax*n+4*n+12*mmax*mmax+12*mmax];
    
    iprint = 1;
    m=5;
    //int i=;
	integer nn=n;
    strncpy(task, "START",5);
    memset(task+5, ' ', sizeof(task)-5);
    while(1) {
	setulb_(&nn,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,csave,lsave,isave,dsave,60,60);
	 if(task[0]=='F') {
               /*FG*/
	     f = func.evaluateFG(x,g,n);
	     //for(int j=0; j<n; j++)
		 //cout << x[j]<<";   ";
	     //cout << endl;
	     continue;
	 }
	 if(task[0]=='N')     /*NEW_X*/ 
	     if (isave[30] > maxit) {
			 cout << "Warning: Maximum iteration limit reached in BFGS \n";
		 break;
	     }
         else continue;                                                       
	 task[59] = 0;
	 cout << "BFGS terminated with message: " << task << endl;
	 break;
    }
    //delete wa;
    delete []iwa;
    delete []g;
    delete []wa;
    return f;
}



