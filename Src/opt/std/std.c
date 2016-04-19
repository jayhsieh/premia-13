#include "std.h"
#include "error_msg.h"

extern Option OPT(CallEuro);
extern Option OPT(CallSpreadAmer);
extern Option OPT(CallSpreadEuro);
extern Option OPT(DigitAmer);
extern Option OPT(DigitEuro);
extern Option OPT(PutAmer);
extern Option OPT(PutEuro);
extern Option OPT(CallAmer);
Option* OPT(family)[]=
{
  &OPT(CallEuro),
  &OPT(PutEuro),
  &OPT(CallSpreadEuro),
  &OPT(DigitEuro),
  &OPT(CallAmer),
  &OPT(PutAmer),
  &OPT(CallSpreadAmer),
  &OPT(DigitAmer),
  NULL
};
