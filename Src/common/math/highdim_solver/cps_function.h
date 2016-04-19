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
#ifndef FUNCTION_H
#define FUNCTION_H

#include "cps_types.h"

#define MAX_PARAMETERS 15

struct function_t {

	const void	*args;
	double 	parameters[MAX_PARAMETERS];
	double	(*body)(const function *, const grid_node *);
};


double cps_function_evaluate(const function *, const grid_node *);
#endif

