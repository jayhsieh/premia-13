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
#ifndef PDE_TERM_H
#define PDE_TERM_H

#include "cps_types.h"

#define UX_TERM 	0xD1
#define UY_TERM 	0xD2
#define UXX_TERM 	0xD3
#define UYY_TERM	0xD4
#define UXY_TERM	0xD5
#define U_TERM		0xD6

struct pde_term_t {

  const function 		*function_factor;
	unsigned int			order;
	unsigned int 			type;

  stencil_operator 	*st_operator;
	stencil 					*generated_stencil;
};

int pde_term_create(pde_term **,int, function *, stencil_operator *);
int pde_term_destroy(pde_term **); 
int pde_term_set_function_factor(pde_term *, const function *);
int	pde_term_set_stencil_operator(pde_term *, stencil_operator *);
int pde_term_create_stencil(pde_term *, const grid *);
#endif

