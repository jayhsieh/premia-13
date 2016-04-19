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
#ifndef PDE_H
#define PDE_H

#include "cps_types.h"

#define PDE_MAX_TERMS 10

struct pde_t {

	pde_term	*terms[PDE_MAX_TERMS];
	const function	*source_term;
	pde_integral_term	*integral_term;
	unsigned int terms_count;
	unsigned int item;
};

int pde_create(pde **);
int pde_destroy(pde **);
int pde_add_term(pde *, pde_term *);
int pde_set_source_term(pde *, const function *);
int pde_set_integral_term(pde *, pde_integral_term *);
int pde_has_source_term(const pde *);
int pde_has_integral_term(const pde *);

int pde_term_start(pde *);
int pde_term_forth(pde *);
int pde_term_item(pde *, pde_term **);
int pde_term_after(pde *);

int pde_next_term(pde *, pde_term **);
#endif

