/***************************************************************************
                          pricer_Avellaneda.h  -  description
                             -------------------
    begin                : Thu Jan 15 2004
    copyright            : (C) 2004 by messaoud
    email                : marouen.messaoud@inria.fr
 ***************************************************************************/

#ifndef PRICER_H
#define PRICER_H

#include "tree.h"
#include "inout.h"

void loadTree(char *name_in_tree , double *pt_S_0, double *pt_r, double *pt_q , int *pt_N ,double *pt_sigma_0,
double*pt_sigma_min, double *pt_sigma_max, double *pt_sigma_bar, double *pt_T_max,char **pt_name_in_data);


double optimPrice(double K, double T, Parametre param , char *OptionType);

#endif