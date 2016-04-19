/***************************************************************************
                          rollback.h  -  description
                             -------------------
    begin                : ven mai 2 2003
    copyright            : (C) 2003 by marouen
    email                : marouen.messaoud@inria.fr
 ***************************************************************************/

#ifndef ROLLBACK
#define ROLLBACK


#include "tree.h"

int define_source( double *source , Vect_option option , Trinomial_tree arbre ,\
                              Vecteur indice , double *lamda , double r ) ;


void define_G_source(double *Gsource , double Ki , char *payoff  , Trinomial_tree arbre , int rank  );


double PHI( double x , /*sigma prior*/double prior , /*sigma min*/ double minp , /*sigma max*/ double maxp ) ;

double GPHI( double x , /*sigma prior*/double prior , /*sigma min*/ double minp , /*sigma max*/ double maxp ) ;

double costFunction ( double *lamda ) ;

void gradCostFunction ( double *lamda , double *grad ) ;

Trinomial_tree getPriorTree(/*vecteur des differentes maturites ordonnees */ Vecteur ,Vecteur *indice , Parametre param);

double price( double K , double S, double T , Parametre param);

double optimPrice(double K, double T, Parametre param , char *OptionType);

#endif


