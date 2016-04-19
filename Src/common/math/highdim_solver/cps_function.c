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
#include <stdlib.h>
#include "cps_utils.h"
#include "cps_function.h"
#include "cps_assertions.h"


double cps_function_evaluate(const function *f, const grid_node *n){
  /* evaluate function in given node */
    double value = f->body(f,n);
	REQUIRE("func_not_null",(f != NULL));
	return value;
}

