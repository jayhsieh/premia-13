#ifndef  _PAD_H
#define _PAD_H

#include  "optype.h"
#include  "var.h"
#include  "chk.h"
#include  "numfunc.h"
#include "option.h"

#define TYPEOPT PAD

/*PathDep Option*/
typedef struct TYPEOPT{
  VAR                                Maturity;                                  
  VAR      PayOff;	/*		The Payoff is phi(stock,path_dep)	*/
  VAR		 PathDep;	/* The PathDep functional definition: 

				new_path-dep=psi(PathDep->Par,stock,time)

				where:

				starting_date is in PathDep->Par[0],
				final_date is in PathDep->Par[1],
				frequency  is in PathDep->Par[2],
				initial_path_dep is in PathDep->Par[3],
				current_path_dep is in PathDep->Par[4]
								
				!!!!!WARNING!!!!!
				Wether the pathdep is backard/forward
				should be tested in ChkOpt
				*/

  VAR					  MinOrElse;	/* cf supra*/   
  VAR					  EuOrAm;
  VAR					  PartOrTot;	/* Partial or total pathdep:

							a partial pathdep is specified
							by starting_date, final_date*/

  VAR					  ContOrDisc;	/*Continuous or Discrete:
							  a discrete pathdep is specified 
                                         by frequency (regular sampling) */
  /* /\*Cliquet options*\/
   * VAR Fg;
   * VAR Cg;
   * VAR Fl;
   *  VAR Cl; */

} TYPEOPT;


/*MinOrElse*/
#define MINIMUM 0
#define MAXIMUM 1
#define AVERAGE 2
 
int OPT(Get)(int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(FGet)(char **InputFile,int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(Show)(int user,Planning *pt_plan,Option *opt, Model *mod);      
int OPT(Check)(int user,Planning *pt_plan,Option *opt);    

#endif
