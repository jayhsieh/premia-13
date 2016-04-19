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
#include <stdio.h>
#include <stdarg.h> 
#include "cps_debug.h"

#ifdef DEBUG
void print_debug(int level, char *fmt, ... )
{
    va_list args;
    char message[MAX_DEBUG_STRING];

    va_start(args, fmt);
    vsnprintf(message, MAX_DEBUG_STRING, fmt, args);
    va_end(args);
		fprintf(stderr,"DEBUG");
		switch(level){
			case DEBUG_INFO:	
				fprintf(stderr,"(info)");
				break;
			case DEBUG_WARN:	
				fprintf(stderr,"(warn)");
				break;
			case DEBUG_ERROR:	
				fprintf(stderr,"(error)");
				break;
		}
    fprintf (stderr, ": %s\n",message);
}
#else
void print_debug(int level, char *fmt, ...) {} ;
#endif

