#ifndef _PREMIA_API_IMPORT_H_INCLUDED_
#define _PREMIA_API_IMPORT_H_INCLUDED_

/*! \file import.h
	\brief Contains imported from premia declarations
*/

//#include <direct.h>

extern "C" {
	
	#include "optype.h"
  #include "enums.h"
  #include "var.h"
	#include "premia_obj.h"
	int InitVar();
	extern int *true_typeV;
	extern int g_dup_printf;
	extern FILE * out_stream;
	extern FILE * g_dup_file;
	extern Model  BSND_model;
    extern Model  COPULA_model;
	extern Family STDND_family;
}

#endif
