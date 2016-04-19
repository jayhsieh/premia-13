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
#ifndef _DEBUG_H
#define _DEBUG_H

#include <stdio.h>
#include <stdarg.h> 

#define MAX_DEBUG_STRING 512

#define DEBUG_INFO 	0xD1
#define DEBUG_WARN 	0xD2
#define DEBUG_ERROR 0xD4

#ifdef CPS_DEBUG
  #ifdef __GNU__
    #define PRINT_DEBUG(f, s...) print_debug(f, ## s)
  #else
        #ifdef WIN32
                #define PRINT_DEBUG(...) { print_debug(__VA_ARGS__); }
        #else
                #define PRINT_DEBUG(f, ...) print_debug(f, __VA_ARGS__)
        #endif
  #endif
        #define DUMP_DEBUG(b, l) dump_debug(b, l)
#else
  #ifdef __GNU__
    #define PRINT_DEBUG(f, s...)
  #else
        #ifdef WIN32
                #define PRINT_DEBUG(...)
        #else
                #define PRINT_DEBUG(f, ...)
        #endif
  #endif
#endif

void print_debug(int, char *fmt, ... );
#endif

