#include "stdnd.h"

extern Option OPT(CallBasketEuro_nd);
extern Option OPT(PutBasketEuro_nd);
extern Option OPT(CallBasketAmer_nd);
extern Option OPT(PutBasketAmer_nd);
extern Option OPT(CallMaximumAmer_nd);
extern Option OPT(PutMinimumAmer_nd);
extern Option OPT(CallMultiSpreadEuro_nd); 
extern Option OPT(PutMultiSpreadEuro_nd); 

extern Option OPT(CallGeomAmer_nd);
extern Option OPT(PutGeomAmer_nd);

Option* OPT(family)[]={
  &OPT(CallBasketEuro_nd),
  &OPT(PutBasketEuro_nd),
  &OPT(CallBasketAmer_nd),
  &OPT(PutBasketAmer_nd),
  &OPT(CallMaximumAmer_nd),
  &OPT(PutMinimumAmer_nd),
  &OPT(CallMultiSpreadEuro_nd),
  &OPT(PutMultiSpreadEuro_nd),
  &OPT(CallGeomAmer_nd),
  &OPT(PutGeomAmer_nd),
  NULL
};
