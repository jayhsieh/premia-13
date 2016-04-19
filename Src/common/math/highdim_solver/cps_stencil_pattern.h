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
#ifndef STENCIL_PATTERN_H
#define STENCIL_PATTERN_H

#include "cps_types.h"

#define GLOC_BOUNDARY 0x1a
#define GLOC_INTERNAL 0x1b
#define GLOC_EXTERNAL 0x1c

struct stencil_application_t {

	double value;
	unsigned int position;
	unsigned int order;
	unsigned short grid_location;
};

struct stencil_pattern_t {

	unsigned int count;
	unsigned int cursor;
	stencil_application *application[MAX_STENCIL_SIZE];
};

int stencil_pattern_create(stencil_pattern **);
int stencil_pattern_destroy(stencil_pattern **);
int stencil_pattern_put(stencil_pattern *, unsigned int, stencil_application *);
/* iterators */
int stencil_pattern_start(stencil_pattern *);
int stencil_pattern_after(const stencil_pattern *);
int stencil_pattern_forth(stencil_pattern *);
int stencil_pattern_item(const stencil_pattern *, stencil_application **);
/* stencil_application */
int stencil_application_create(stencil_application **);
int stencil_application_destroy(stencil_application **);
int stencil_application_is_boundary(const stencil_application *);
int stencil_application_is_external(const stencil_application *);
int stencil_application_is_internal(const stencil_application *);
int stencil_application_set_boundary(stencil_application *);
int stencil_application_set_external(stencil_application *);
int stencil_application_set_internal(stencil_application *);
int stencil_application_set_order(stencil_application *, unsigned int);
#endif

