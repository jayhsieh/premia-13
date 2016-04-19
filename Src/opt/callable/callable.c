#include "callable.h"
#include "error_msg.h"

extern Option OPT(NoCall);
extern Option OPT(StdCall);
extern Option OPT(PathDepCall);
extern Option OPT(HighlyPathDepCall);
extern Option OPT(IntermittentCall);

Option* OPT(family)[]=
{
  &OPT(NoCall),
  &OPT(StdCall),
  &OPT(PathDepCall),
  &OPT(HighlyPathDepCall),
  &OPT(IntermittentCall),
  NULL
};
