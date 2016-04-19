#include "std2dg.h"
#include "error_msg.h"

extern Option OPT(CallSpread2dEuro);
extern Option OPT(PutSpread2dEuro);

Option* OPT(family)[]=
{
  &OPT(CallSpread2dEuro),
  &OPT(PutSpread2dEuro),
  NULL
};
