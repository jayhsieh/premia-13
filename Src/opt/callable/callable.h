#ifndef  _STD_H
#define _STD_H

#include  "optype.h"
#include  "var.h"
#include  "option.h"
#include  "chk.h"

#define TYPEOPT CALLABLE

typedef struct TYPEOPT
{         
  VAR Strike; /* Final Strike, if no anticipative exercise */
  VAR Maturity; /* Maturity */
  VAR Coupon; /* Nominal coupon rate */
  VAR Recovery; /* Nominal Recovery */
  VAR PutStrike; /* Strike if the holder exercises */
  VAR CallStrike; /* Strike if the seller exercises */
  VAR LowerBarrier; /* Lower Barrier */
  VAR UpperBarrier; /* Upper Barrier */
  VAR Window; /* Window for Parisian type condition */
  VAR Period; /* Moving period over which we check the condition */
} TYPEOPT;


#endif
