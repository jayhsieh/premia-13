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
#ifndef STENCIL_OPERATOR_H
#define STENCIL_OPERATOR_H

#include "cps_types.h"

#define STENCIL_OP_UXX 	0xA7
#define STENCIL_OP_UYY 	0xA6
#define STENCIL_OP_UXY 	0xA5
#define STENCIL_OP_UX		0xA4
#define STENCIL_OP_UY		0xA3
#define STENCIL_OP_U		0xA2

struct stencil_operator_t {

	unsigned int type;
	unsigned int is_applied;
	
	stencil* 	applied_stencil;
	stencil* 	(*apply)(const pde_term *, const grid *);
};

int stencil_operator_create(stencil_operator **, int);
int stencil_operator_destroy(stencil_operator **);
int stencil_operator_apply(stencil_operator *, const pde_term *, const grid *);
#endif

