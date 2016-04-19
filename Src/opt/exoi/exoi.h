#ifndef  _EXOi_H
#define _EXOi_H

#include  "optype.h"
#include  "var.h"
#include  "chk.h"
#include  "numfunc.h"
#include "option.h"

#define TYPEOPT EXOi

typedef struct TYPEOPT
{
    VAR PayOff;
    VAR FirstExerciseDate;
    VAR LastPaymentDate;
    VAR ResetPeriod;
    VAR Nominal;
    VAR Spread;
    VAR Cap;
    VAR Strike;
    VAR Gearing;
    VAR Floor;
    VAR FixedRate;
    VAR LowerRangeBound;
    VAR UpperRangeBound;
    VAR CMSMat1;
    VAR CMSMat2;
} TYPEOPT;

/*variables that define the coupon :*/
// CallableCappedFloater=(Spread, Cap)
// CallableInverseFloater=(Cap, Strike, Gearing, Floor)
// CallableRangeAccrual=(FixedRate, LowerRangeBound, UpperRangeBound)


int OPT(Get)(int user,Planning *pt_plan,Option *opt, Model *mod);
int OPT(FGet)(char **InputFile,int user,Planning *pt_plan,Option *opt, Model *mod);
int OPT(Show)(int user,Planning *pt_plan,Option *opt, Model *mod);
int OPT(Check)(int user,Planning *pt_plan,Option *opt);

#endif
