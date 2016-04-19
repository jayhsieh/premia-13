#include "stdf.h"

extern Option OPT(InflationIndexedCap);
extern Option OPT(InflationIndexedCaplet);
extern Option OPT(YearOnYearInflationIndexedSwap);

Option* OPT(family)[]=
{
  
  &OPT(InflationIndexedCap),
  &OPT(InflationIndexedCaplet),
  &OPT(YearOnYearInflationIndexedSwap),
  NULL
};
