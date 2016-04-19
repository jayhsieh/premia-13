#include "stda.h"

extern Option OPT(EquityLinkedSurrenderEndowment);
extern Option OPT(Gap);

Option* OPT(family)[]=
{
  &OPT(EquityLinkedSurrenderEndowment),
  &OPT(Gap),
  NULL
};
