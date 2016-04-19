#ifndef LMM_NUMERICAL_H 
#define LMM_NUMERICAL_H


#include"lmm_header.h"
#include"lmm_volatility.h"
#include"lmm_libor.h"
#include"lmm_random_generator.h"
#include"lmm_zero_bond.h"
#include"lmm_mathtools.h"


extern double maxi(double a,double b);
extern double ppos(double x);
extern int Set_to_Zero(double *x,int dim);
extern int evolutionUnderSpotMeasure(RandomGenerator *ptRand, Libor* ptLibOld , Libor* ptLibNew, Volatility *ptVol, double dt , double t);
extern int evolutionUnderForwardMeasure(RandomGenerator *ptRand, Libor* ptLibOld , Libor* ptLibNew, Volatility *ptVol, double dt , double t);
extern void computeNumeraire(char*nameOfthemeasure,Libor* ptLib,Swaption* ptSwpt,double *Numeraire,int j,int k,double auxspot);
extern void Name_To_Measure(char *ErrorMessage, char *name,
		     int (**computeEvol)(RandomGenerator* ptRand,Libor* ptLibOld,Libor* ptLibNew,Volatility* ptVol,double dt,double t));
extern int computeSwaptionPriceSpot( int numberMonteCarloSim, int numberTimeStep , double dt ,Libor *ptLib , Swaption *ptSwpt, Volatility *ptVol );
extern int computeSwaptionPriceForward( int numberMonteCarloSim, int numberTimeStep , double dt ,Libor *ptLib , Swaption *ptSwpt, Volatility *ptVol );
extern int computeZeroBond( int numberMonteCarloSim, int numberTimeStep , double dt ,Libor *ptLib , Volatility *ptVol , ZeroBond * ptZb);
extern int computeTrajectories( int numberTimeStep , double dt ,HistLibor *ptHistLib ,  Volatility *ptVol );

////

extern double ps(int n, double *u, double *v); /*scalar product of 2 vectors*/
extern double BSFormula(Swaption *ptSwpt, Libor* ptLib, double evalTime, double blackVol);
extern double maximum(double a, double b);


#endif
