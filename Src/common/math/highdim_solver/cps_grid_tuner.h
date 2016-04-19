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
#ifndef GRID_TUNER_H
#define GRID_TUNER_H

#include "cps_types.h"
#define MAX_TUNERS 4

#define GENERIC_TUNER 	0
#define EXPLICIT_TUNER 	1
#define IMPLICIT_TUNER	2
#define RESCALE_TUNER 	3

typedef int (*grid_tuner_proc)(grid_tuner *, grid *);

struct grid_tuner_t {

	void *argument;
	grid_tuner_proc tuners[MAX_TUNERS];	
};

int grid_tuner_create(grid_tuner **);
int grid_tuner_destroy(grid_tuner **);
int grid_tuner_set_tuner(grid_tuner *, int, grid_tuner_proc);
int grid_tuner_set_argument(grid_tuner *, void *);
int grid_tuner_apply(grid_tuner *, int, grid *);
#endif

