#include "stdg.h"
#include "error_msg.h"

extern Option OPT(CallFuture);
extern Option OPT(PutFuture);
extern Option OPT(Swing);

Option* OPT(family)[]=
{
  &OPT(CallFuture),
  &OPT(PutFuture),
  &OPT(Swing),
  NULL
};
