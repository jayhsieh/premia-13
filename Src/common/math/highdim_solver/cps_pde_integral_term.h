/*********************************************************
*   CPS - A simple C PDE solver                          * 
*                                                        *
*   Copyright (c) 2007,                                  * 
*     Maya Briani       <m.briani@iac.rm.cnr.it>,        *        
*     Francesco Ferreri <francesco.ferreri@gmail.com>,   * 
*     Roberto Natalini  <r.natalini@iac.rm.cnr.it>,      *
*     Marco Papi        <m.papi@iac.rm.cnr.it>           *
*                                                        * 										
*********************************************************/
#ifndef PDE_INTEGRAL_TERM_H
#define PDE_INTEGRAL_TERM_H

#include "cps_types.h"
#include "laspack/highdim_vector.h"

struct pde_integral_term_t {

	double lambda;
	double alpha;
	double m;

	const grid *source_grid;
};

int pde_integral_term_create(pde_integral_term **);
int pde_integral_term_destroy(pde_integral_term **); 
int pde_integral_term_set_lambda(pde_integral_term *, double);
int pde_integral_term_set_m(pde_integral_term *, double);
int pde_integral_term_set_alpha(pde_integral_term *, double);
int pde_integral_term_set_grid(pde_integral_term *, const grid *);
double pde_integral_term_evaluate(const pde_integral_term *, const grid_node *, Vector *);
#endif

