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
#ifndef BOUNDARY_DESCRIPTION_H
#define BOUNDARY_DESCRIPTION_H

#include "cps_types.h"
#include "cps_dimensions.h"

#define BOUNDARY_NONE 	0xB0
#define BOUNDARY_LEFT 	0xB1
#define BOUNDARY_RIGHT 	0xB2
#define BOUNDARY_INITIAL 0xB3

struct boundary_description_t {

	function	*initial;
	function	*left[MAX_DIMENSIONS];
	function	*right[MAX_DIMENSIONS];
};

int boundary_description_create(boundary_description **);
int boundary_description_set_left(boundary_description *, unsigned int, function *);
int boundary_description_set_right(boundary_description *, unsigned int, function *);
int boundary_description_set_initial(boundary_description *, function *);
double boundary_description_evaluate(boundary_description *, const grid *, const grid_node *);
int boundary_description_destroy(boundary_description **);
#endif

