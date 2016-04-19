#ifndef LMM_BERMUDAPRICE_H
#define LMM_BERMUDAPRICE_H

#include"lmm_header.h"


/*Main Routine*/
extern int computeBermudeanSwaption(long numberMCPaths,int numberTimeStep,
					Libor* ptLib, Swaption* ptSwpt, Volatility* ptVol,
					int RegrBasisDimension,double PayOff_As_Regressor,
					char* Basis_Choice,char* Measure_Choice,char Explanatory_Var);


/*allocation/liberation routines for local variable necessary to price bermudan swaptions*/

extern void mallocBermudaVar(double **RegrVar,double **Res,double** Brownian,double** SwapPrices,double** Numeraire,double **FP,
			long numberMCPaths, int RegrVarDimension,int RegrBasisDimension,int numberOfExerciseDates,int Brown_factors,
			     char* ErrorMessage);
extern void freeBermudaVar(double **RegrVar,double **Res,double** Brownian,double** SwapPrices,double** Numeraire);

extern void Regression(long NumberMCPaths,int numberOfExerciseDates,
		       int RegrBasis_Dimension, 
		       int X_Dimension, 
		       int Swap_Entry_Time,
		       int PayOff_As_Regressor,
		       double* X,double* FP,double* Swap_Prices,
		       double* Res,
		       char* ErrorMessage );

extern void initStateVector(double *X,double *Brownian,double* SwapPrices,double *Numeraire,int j,long numberMCpaths,int numberOfExerciseDates,int numberOfFactors,char Rflag);
 
extern double lmm_swaption_payer_bermudan_LS_pricer(float tenor , int numberTimeStep, int numFac, double swaptionMat , double swapMat ,  double payoff_as_Regressor ,long  numberMCPaths , int Regr_Basis_Dimension , char* basis_name , char* measure_name , char Explanatory , double strike);


#endif
