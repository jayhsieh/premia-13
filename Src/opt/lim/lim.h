#ifndef  _LIM_H
#define _LIM_H

#include  "optype.h"
#include  "var.h"
#include  "chk.h"
#include  "numfunc.h"
#include "option.h"

#define TYPEOPT LIM

/*Limit Option// Single barrier*/

typedef struct TYPEOPT{   
  /* setable */
  VAR Maturity;   
  VAR Limit;      /*The Limit definition:
		   * starting_date is in Limit->[0],
		   * final_date is in Limit->Par[1],
		   * frequency is in Limit->Par[2],
		   * the value of the Limit in case of a constant limit is in Limit->Par[3]
		   * Parisian delay is in Limit->Par[4],
		   * !!!!!WARNING!!!!!
		   * Wether the limit is backard/forward
		   * should be tested in ChkOpt
		   */
  VAR PayOff;
  VAR Rebate;
  /* non setable */
  VAR OutOrIn;
  VAR Parisian;
  VAR DownOrUp;
  VAR RebOrNo;
  VAR EuOrAm; 
  VAR PartOrTot; /* Partial Or Total limit
		  * a partial limit is specified
		  * by starting_date, final_date
		  */
  VAR ContOrDisc;/*Continuous or Discrete:
		  * a discrete limit is specified 
		  * by frequency (regular sampling) 
		  */ 
  VAR ConstLim;  /*YES for constant limit*/

} TYPEOPT;

int OPT(Get)(int user,Planning *pt_plan,Option *opt, Model *mod); 
int OPT(FGet)(char **InputFile,int user,Planning *pt_plan,Option *opt, Model *mod);   
int OPT(Show)(int user,Planning *pt_plan,Option *opt, Model *mod);
int OPT(Check)(int user,Planning *pt_plan,Option *opt); 

#endif



