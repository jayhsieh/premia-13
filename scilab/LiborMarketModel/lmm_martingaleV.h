#ifndef LMM_MARTINGALEV_H
#define LMM_MARTINGALEV_H


#include"lmm_header.h"
#include"lmm_volatility.h"
#include"lmm_libor.h"
#include"lmm_numerical.h"
#include"lmm_random_generator.h"
#include"lmm_products.h"

/*Structure of the variables V_n which are martingales for the spot measure and whose discretisations stay martingale*/
typedef struct
{
  int NbMaturity;     /* Number of maturities T0, T1, T2 ..., TNbMaturity */
  double tenor;       /* constant tenor between 2 maturities */
  double *maturity;   /* Maturities of all underlaying libor asset*/
  double *Vvalue_0;   /* Initial value of the dicrete martingale variable V*/ 
  double *Vvalue_t;   /* Value after a evolution till time t resulting of the function drawMartingaleV() */
  
}MartingaleV;


extern int mallocMartingleV(MartingaleV **ptV, Libor ptL);
/*Allocate with LIBORS values the Vvalue_0 and Vvalue_t (equal to the Vvalue_0), V has one more value than L*/
extern int initMartingaleV(MartingaleV *ptV);
/*Set the Vvalue_t equal to the Vvalue_0*/
extern int freeMartingaleV(MartingaleV **ptV);
/*Free the memory of the structure MartingaleV*/
extern int printMartingaleV(MartingaleV *ptV);
/*Print the main datas of the object under ptV*/
extern int drawMartingaleV(MartingaleV *ptV, Volatility *ptVol, RandomGenerator *ptW, double dt, double t0,double t);
/*Function which draws the variable V from t0 till time t according to a log Euler schem */
extern int SetLiborWithMartingaleV(Libor *ptL, MartingaleV *ptV);
/*Set the values of the Libors values corresponding to a the values of ptV (according to the formula that links both)*/
extern int SetMartingaleVWithLibor(MartingaleV *ptV, Libor *ptL);
/*Set the values of ptV corresponding to a the value of the Libors (according to the formula that links both)*/
extern double computeSwaptionSpotV(int Nbmontecarlo,double dt,Libor *ptLib,Swaption *ptSwpt,Volatility *ptVol);
/*Return the value of the swaption defined under ptSwpt by creating (thanks to ptLib) a martingaleV object drawn Nbmontecarlo times*/
extern int computeCapletSpotV(int Nbmontecarlo,double dt,Libor *ptLib, double K, Volatility *ptVol, double *Caps);
/*Compute all the caplets prices stored in Caps. double pointer Caps must have the correct allocated size= Nbmaturity (Cap[0] is nothing, Cap[1] maturing at T1 paying at T2, Cap[2]...*/
extern int computeBiasSpotV(int Nbmontecarlo,double dt,Libor *ptLib1,double K, Volatility *ptVol, double *Bias);
/*Compute the caplets bias (with the martingaleV drawing and montecarlo) compared to the BS formula (stored in Bias wich must have to correct allocated size=NbMaturity) Not really correctly tested, DONT USE IT*/
extern int computeZCBondSpotV(int Nbmontecarlo,double dt,Libor *ptLib, double maturity, Volatility *ptVol);
/*Compute the ZCbonds bias (with the martingaleX drawing and montecarlo) compared to exact formula print in file "zcbonderror.dat" Not really correctly tested, DONT USE IT*/
extern int printdoubles(char *filename, int n, double *D);

int lmm_caplet_spotV_pricer(double *caplets_price ,  double* maturities , int nb_mat , int nb_MC , int nb_factors , int nb_time_step , double strike , double tenor);
int lmm_swaption_payer_spotV_pricer(double *swaption_price ,  double swaption_maturity , double swap_maturity    , int nb_MC , int nb_factors , int nb_time_step , double strike , double tenor);

#endif




