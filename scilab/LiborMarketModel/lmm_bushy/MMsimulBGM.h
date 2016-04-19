#ifndef MMSIMULBGM_H
#define MMSIMULBGM_H
#include <memory.h>
//#include <stdio.h>
//#include<stdlib.h>
//#include"math.h"

#include <stdio.h>
#include <stdlib.h>
//#include <ctype.h>
#include <math.h>
#include <assert.h>

#define MAX(x,y) ((x>y) ? (x):(y))
long NSteps,NRates,NBranches,NumeraireIndex;
//LogShiftOfBranch[h][k][i] contains contains -1/2 cii* dt +bik*sqrt(dt) for the time step from th to th+1 , 
//the array element C[h][i][k] holds the associated covariance matrix entry . cik for the time step,

double ***C;

double *** LogShiftOfBranch,***Bbar;
double ** EvolvedFra,**S,** A, **Mattemp;
double * Tau, **mu_dT;
double period,Ts,strike;
double B0;//ZeroCoupon
int flag,bushyprixcalculer;
int type;// 1 pour (S-K)+ et -1 pour (K-S)+
int typeOption; // 1 Swaption ; 0 si (Caplet - Floorlet);
FILE* fic;
int flagBushy,flagRotate;

//parametre pour Bushy tree
long bushytimeStep,bushyBranches;
long *bushyInd;	  // de taille NSteps/bushytimeStep +1;
double **bushySwap;// de taile( NBranches^bushytimeStep)*(NSteps/bushytimeStep +1)
double **bushyprix,bushyprixnoeud;// de taille NBranches^bushytimeStep
double *y,*x;//de taille bushyBranches
/*------------------------ the main function that start the recusion and-----------------------------------*/
/*-------------------------------- that return the non actualized payoff-----------------------------------*/
double Recurse(long h);
double RecurseBushy(long h);
double RecurseCaplet(long h);
double RecurseBushyCaplet(long h);
double Interface();
/*-------------construct the vector, allocate the memeory , -----------------------------------------------*/
/*------------------call "double Recurse(long h)" and return the price ----------------------------------- */
double Price();

/* return the payoff at time h*/
double Intrinsic(long h);
double CheckForEarlyExercise(long h,double esp );
/*-----------------set the LogShiftOfBranch at time h--------------------------------------------*/
int SetLogShiftOfBranch(long h);
/*-----------------set S the shift transition matrix---------------------------------------------*/
int SetS();
int SetBbar(long h);
int setCovarMatrix(double **C ,long h);
/* ---------------------------return the decorrelated matrix in C--------------------------------*/
int DecorrelateCovarianceMatrix(double **M, int NRate,int NFactor, double ** C,int h);
int ConstructCovarianceMatrix(int NRate,int NFactor, double **C,int h);
/*cholesky methods*/
void Cholesky(double **M, int dim);
int MyCholesky(double **M,double **B, int dim);
int AfficheMat(double **M,int diml,int dimc);
int RotateSimplex(int NFactor);
double costFunction ( double *lamda ) ;
void gradCostFunction ( double *lamda , double *grad ) ;
void matriceprod(double ** A,double **B,double **res,long nl, long nc);


/*----------- allocate memory to C the covariance matrix of size (NSteps *(NRates*NRates))  ---------------*/
/* si la matrice de covaraince est constante pas besoin de la dimension NSteps ----------------------------*/
int C_AllocatMem();
/*Allocate a memory to matrix of size[NRates][NBranches-1] will be used to stock the cholesky matrix-------*/
/*----------------------------------------- such that A*A'=C-----------------------------------------------*/
int A_AllocatMem();
/*----------- allocate memory to S the shift transition matrix of size (NBranches*(NFactor=NBranches-1))---*/
/* ------------ contains the up down midle.... of the tree ------------------------------------------------*/
int S_AllocatMem();

/*------------------Tau[] de taille[NSteps+1]; c une approximation grossaire ---------------------*/
/* ----------normalement pour les abres il faut rafiner le pas de temps-------------------------*/
int Tau_AllocatMem();



/*------------------ de taille[NRates];  contient les drifts des Libor --------------------------*/

int mu_dT_AllocatMem();

/*------------------ de taille[NSteps+1][NRates];  contient les Libor dans le temps--------------------------*/
int EvolvedFra_AllocatMem();
/*------------------ de taille[NSteps][NBranches];  contient les gaussiennes*variance dans le temps--------*/
/*------------------------------Bbar=A*S' -----------------------------------------------------------------*/
int Bbar_AllocatMem();

/*----------- allocate memory to LogShiftOfBranch (NSteps *(NBranches*NRates))  ---------------*/
/* ----- contains at h ,k,i = -0.5*C[i][i]+Bbar[i][k] -----------------------------------------*/
int LogShiftOfBranch_AllocatMem();

int C_FreeMem();
int S_FreeMem();
int A_FreeMem();
int Tau_FreeMem();
int mu_dT_FreeMem();
int	EvolvedFra_FreeMem();
int	Bbar_FreeMem();
int	LogShiftOfBranch_FreeMem();

int Allocate_All_memory();
int Free_All_memory();
#endif
