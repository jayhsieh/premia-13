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
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "cps_assertions.h"

void assertion_callback(assert_type type, const char *file, int line, const char *func, const char *tag, const char *cond, const char *addinfo){ 

	const char *str_type;

  switch (type) {
  case assert_require:
    str_type = "PRECONDITION";
    break;
  case assert_ensure:
    str_type = "POSTCONDITION";
    break;
  case assert_check:
    str_type = "CHECK";
    break;
  default:
    str_type = "UNKNOWN";
  }

	fprintf(stderr,"\n*** assertion violation (%s) in %s() (%s,%d)\n\n", str_type, func, file, line);
	fprintf(stderr,"%s: (%s)\n", tag, cond); 
	if(addinfo){
		fprintf(stderr," ,%s\n",addinfo);
	}
	else{
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"***\n");
  abort();
}

