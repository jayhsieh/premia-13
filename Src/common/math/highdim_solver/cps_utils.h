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
#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <string.h>
#include "cps_types.h"
#include "cps_assertions.h"
#include "cps_debug.h"

/* utility macros */

#ifdef WITH_GC
#include <gc.h>
#define APP_INIT GC_INIT()
#define MEM_ALLOC GC_MALLOC
#define MEM_FREE  GC_FREE
#else
#define APP_INIT
#define MEM_ALLOC malloc
#define MEM_FREE free
#endif


#define STANDARD_CREATE(PTR,TYPE)							\
	*PTR = (TYPE*)MEM_ALLOC(sizeof(TYPE));			\
	CHECK("ptr_not_null",((*(PTR)) != NULL));				\
  memset((*(PTR)),0,sizeof(TYPE));								

#define STANDARD_DESTROY(PTR)									\
	if((*(PTR)) != NULL)					\
		MEM_FREE((*(PTR)));
#endif

