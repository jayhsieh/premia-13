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
#ifndef PROBLEM_SOLVER_H
#define PROBLEM_SOLVER_H

#include "laspack/qmatrix.h"
#include "laspack/highdim_vector.h"
#include "laspack/itersolv.h"
#include "laspack/operats.h"
#include "laspack/errhandl.h"
#include "laspack/rtc.h"

#include "cps_types.h"

#define SOLVER_MODE_IMP 0xF1
#define SOLVER_MODE_EXP 0xF0

#define SOLVER_ALG_CG			0xA1
#define SOLVER_ALG_GMRES	0xA2
#define SOLVER_ALG_BICGS	0xA3

#define MAX_MAIN_SOLVER_ITERATIONS 20
#define MAX_BACKUP_SOLVER_ITERATIONS 100

#define FULL_CORRECTION 0xC1
#define FAST_CORRECTION 0xC2

struct problem_solver_t {

	int mode;
	int step;
	int algorithm;
	int correction_mode;
	pde_problem *problem;
	QMatrix	Dc,Dn;
	Vector	uc,un,bc;
	IterProcType	iterative_solver;	
};

int problem_solver_create(problem_solver **);
int problem_solver_destroy(problem_solver **);
int problem_solver_setup(problem_solver *, pde_problem *);
int problem_solver_reset(problem_solver *);
int problem_solver_set_mode(problem_solver *, int);
int problem_solver_set_correction_mode(problem_solver *, int);
int problem_solver_set_algorithm(problem_solver *, int);
int problem_solver_step(problem_solver *);
int problem_solver_get_solution_element(problem_solver *, unsigned int, double *);
#endif

