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
#ifndef TYPES_H
#define TYPES_H



typedef struct pde_problem_t pde_problem;
typedef struct problem_solver_t problem_solver;
typedef struct pde_term_t pde_term;
typedef struct pde_integral_term_t pde_integral_term;
typedef struct pde_t pde;
typedef struct grid_t grid;
typedef struct grid_tuner_t grid_tuner;
typedef struct grid_node_t grid_node;
typedef struct function_t function;
typedef struct stencil_t stencil;
typedef struct stencil_pattern_t stencil_pattern;
typedef struct stencil_application_t stencil_application;
typedef struct stencil_operator_t stencil_operator;
typedef struct boundary_description_t boundary_description;
#endif

