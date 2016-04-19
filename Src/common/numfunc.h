#ifndef  _NUMFUNC_H
#define _NUMFUNC_H

double Call(VAR*,double);
double Put(VAR*,double);
double CallSpread(VAR*,double);
double Digit(VAR*,double);
double Zero(VAR*,double);
double Const(VAR*,double);
double ConstLim(VAR*,double);
double DigitSpecialPayoff(VAR *,double spot);

double BestOf(VAR*param,double spot1,double spot2);
double CallMax(VAR*param,double spot1,double spot2);
double Geom(VAR*param,double spot1,double spot2);
double Arim(VAR*param,double spot1,double spot2);
double PutMin(VAR*param,double spot1,double spot2);
double Exchange(VAR*param,double spot1,double spot2);
double Zero2d(VAR*param,double spot1,double spot2);
double Const2d(VAR*param,double spot1,double spot2);
double Call_2arg(VAR*param,double spot1,double spot2);
double Put_2arg(VAR *param,double spot1,double spot2);
double Call_OverSpot2(VAR *param,double spot1,double spot2);
double Put_OverSpot2(VAR *param,double spot1,double spot2);
double Call_StrikeSpot2(VAR*param,double spot1,double spot2);
double Put_StrikeSpot2(VAR*param,double spot1,double spot2);
double CallSpread2d(VAR*param,double spot1,double spot2);
double PutSpread2d(VAR*param,double spot1,double spot2);

double Minimum(VAR*param,double ,double);
double Maximum(VAR*param,double ,double);
double Asian(VAR*param,double,double);


double PutBasket_nd(VAR *param,PnlVect *VStock);
double CallBasket_nd(VAR *param,PnlVect *VStock);
double CallMultiSpread_nd(VAR *param,PnlVect *VStock);
double PutMultiSpread_nd(VAR *param,PnlVect *VStock);
double CallMax_nd(VAR *param,PnlVect *VStock);
double PutMin_nd(VAR *param,PnlVect *VStock);
double PutGeom_nd(VAR *param,PnlVect *VStock);
double CallGeom_nd(VAR *param,PnlVect *VStock);


#endif
