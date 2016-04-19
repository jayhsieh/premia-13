#include <stdio.h>
#include <math.h>

#include "QuasiNewton.h"
#include "BFGSupdate.h"
#include "cholesky.h"
#include "lineSearch.h"
#include "stopping.h"
#include "testing.h"


int main()
{
 /*  testMerite(); */
/*   testWolfeRec_dicho(); */
/*   testWolfeRec_cubic(); */
/*   testcubicMin(); */
  //testWolfe_cubic();
  //testChol();
  //testBFGS();
  //testQuasiNewton(Rosenbrock100,3);
  //testQuasiNewton(WoodProblem,3);
/*   fillBigGaussSol(); */
/*   fillBigGaussStart(); */
/*   testQuasiNewton(BigGaussProblem,3); */
  testQuasiNewton(GaussProblem,1);
  //testQuasiNewton(Quad2,

  return 0;

}
