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
#ifndef STENCIL_H
#define STENCIL_H

#include "cps_function.h"
#include "cps_grid.h"
#include "cps_grid_node.h"

#define MAX_STENCIL_SIZE 9

#define MAX_MODES 2
#define MAX_TIMES 2

#define MODE_EXP 0
#define MODE_IMP 1
#define TIME_CUR 0	
#define TIME_NXT 1

#define XY			0		/* i,j      */
#define	XPY			1		/* i+1,j    */
#define XPYM		2		/* i+1,j-1	*/
#define XYM			3		/* i,j-1		*/
#define XMYM		4		/* i-1,j-1	*/
#define	XMY			5		/* i-1,j		*/
#define XMYP		6		/* i-1,j+1	*/
#define XYP			7		/* i,j+1		*/	
#define	XPYP		8		/* i+1,j+1	*/

struct stencil_t {

	double				weight[MAX_TIMES][MAX_MODES];
	double				factor;
	const 				function *function_factor;
	double				value[MAX_STENCIL_SIZE];
};

int stencil_create(stencil **);
int stencil_destroy(stencil **);
int stencil_set_factor(stencil *, double);
int stencil_set_function_factor(stencil *, const function *);
int stencil_set_value(stencil *, int, double);
int stencil_set_weight(stencil *, int, int, double);
int stencil_apply(stencil *, const grid *, int, int, const grid_node *, stencil_pattern **);
int stencil_evaluate(stencil *, int, int, int, const grid_node *, double *);
#endif

