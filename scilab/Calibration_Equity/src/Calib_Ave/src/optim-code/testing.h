//
// MATHFI Project, Inria Rocquencourt.
// Vincent Barette, June 2002.
//

/* 
** FILE
** *****************************************************************
** 
** testing.h
** 
** 
** PURPOSE
** 
** Some test-problems in order to test optimizer.
** See:
**      - Schittkowski, "More Test Examples for Nonlinear Programming Codes"
** (Springer-Verlag, 1987).
**      - J. E. Dennis and R.B. Schnabel, "Numericals Methods for Unconstrained
** Optimization and Nonlinear Equations" (Prentice-Hall, 1983)
** 
*****************************************************************
*/

#ifndef TESTING
#define TESTING


/* 
** STRUCTURE
** *****************************************************************
** 
** testProblem
**
** DESCRIPTION
** 
** contains:
**    - the size of the problem
**    - the starting point
**    - the cost-function F
**    - the gradient of F
**    - the solution of the problem
**  
*****************************************************************
*/

typedef struct 
{
  int size;
  double* startingPoint;
  double (*costFunction)(double*);
  void (*gradCostFunction)(double*, double*);
  double* solution;
} testProblem;



/* 
** TEST-PROBLEMS
** *****************************************************************
** 
** some testProblems
**
*****************************************************************
*/

extern const testProblem Quadratic;
extern const testProblem Quadratic2;
extern const testProblem Rosenbrock;
extern const testProblem Rosenbrock100;
extern const testProblem WoodProblem;
extern const testProblem GaussProblem;
extern const testProblem BigGaussProblem;


#endif
