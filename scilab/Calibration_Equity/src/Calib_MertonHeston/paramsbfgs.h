#ifndef _PARMSBFGS_
#define _PARMSBFGS_

////////////////////////////////////////////////////////
// Parametres pour le L-BFGS : le BFGS a memoire limitee
// avec bornes
// number of corrections used in the limited memory matrix
#define mg0  5;
// ~ tolerance sur la variation relative (* presicion machine )de f
//: e12 faible precision, e7 moyenne, et e1 tres grande
#define factr0  1.e7;
// tolerace sur la norme du gradient
#define pgtol0  1.e-12  ;
// pour le type d'affichage a chaque iteration
#define iprint0  0;
// nombre maximum d'iterations
#define maxiter0  100;


////////////////////////////////////////////////////////
// Parametres pour le calcul du gradient par differences finies
// l'increment relatif pour FD : x -> x + x*h
#define h0 1.e-6

#endif
