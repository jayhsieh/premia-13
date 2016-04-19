#ifndef LMM_MARTINGALEX_H
#define LMM_MARTINGALEX_H


#include"lmm_header.h"
#include"lmm_volatility.h"
#include"lmm_libor.h"
#include"lmm_numerical.h"
#include"lmm_random_generator.h"
#include"lmm_products.h"

/*Structure of the variables X_n which are martingales for the terminal measure and whose discretisations stay martingale */
typedef struct
{
  int NbMaturity;     /*Number of maturities T_0, T_1, T_2 ..., T_NbMaturity */
  double tenor;       /*Constant tenor between 2 maturities */
  double *maturity;   /*Maturities of all underlaying libor asset*/
  double *Xvalue_0;   /*Initial value of the dicrete martingale variable X (there are NbMatutity value of Xvalue_0)*/ 
  double *Xvalue_t;   /*Value after a diffusion resulting of drawMartingaleX() (there are NbMatutity value of Xvalue_0)*/
  
}MartingaleX;



extern int mallocMartingleX(MartingaleX **ptX, Libor ptL);
/*Allocate with LIBORS values the Xvalue_0 and Xvalue_t (equal to the Xvalue_0), X and L have the same number of value*/
extern int initMartingaleX(MartingaleX *ptX);
/*Set the Xvalue_t equal to the Xvalue_0*/
extern int freeMartingaleX(MartingaleX **ptX);
/*Free the memory of the structure MartingaleX*/
extern int printMartingaleX(MartingaleX *ptX);
/*Print the main datas of the object under ptX*/
extern int drawMartingaleX(MartingaleX *ptX, Volatility *ptV, RandomGenerator *ptW, double dt, double t0,double t);
/*Function which draws the variable X from time t0 till time t according to a log Euler schem */
extern int SetLiborWithMartingaleX(Libor *ptL, MartingaleX *ptX);
/*Set the values of the Libors corresponding to a the value of ptX (according to the formula that links both*/
extern double computeSwaptionTerminalX(int Nbmontecarlo,double dt,Libor *ptLib,Swaption *ptSwpt,Volatility *ptVol);
/*Return the value of the swaption defined under ptSwpt by creating (thanks to ptLib) a martingaleX object and drawind it Nbmontecarlo times*/
extern int computeCapletTerminalX(int Nbmontecarlo,double dt,Libor *ptLib, double K, Volatility *ptVol, double *Caps);
/*Compute all the caplets prices stored in Caps. double pointer Caps must have the correct allocated size= Nbmaturity (Cap[0] is nothing, Cap[1] maturing at T1 paying at T2, Cap[2]...*/
extern int computeBiasTerminalX(int Nbmontecarlo,double dt,Libor *ptLib1,double K, Volatility *ptVol, double *Bias);
/*compute the caplets bias (with the martingaleX drawing and montecarlo) compared to the BS formula (stored in Bias wich must have to correct allocated size=NbMaturity)*/
extern int computeZCBondTerminalX(int Nbmontecarlo,double dt,Libor *ptLib, double maturity, Volatility *ptVol);
/*compute the ZCbonds bias (with the martingaleX drawing and montecarlo) compared to exact formula print in file "zcbonderror.dat"*/
extern int printdoubles(char *filename, int n, double *D);


int lmm_caplet_terminalX_pricer(double *caplets_price ,  double* maturities , int nb_mat , int nb_MC , int nb_factors , int nb_time_step , double strike , double tenor);
int lmm_swaption_payer_terminalX_pricer(double *swaption_price ,  double swaption_maturity , double swap_maturity   ,  int nb_MC , int nb_factors , int nb_time_step , double strike , double tenor);

#endif




