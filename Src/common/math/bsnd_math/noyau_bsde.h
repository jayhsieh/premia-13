
#ifndef _NOYAU_BSDE_H
#define _NOYAU_BSDE_H

#include "pnl/pnl_matrix.h"

#define EPS_BSDE 0.0000001

extern void init_alpha(int d);
extern double operateur_noyau(PnlMat *X, PnlVect *T_alea, double s,
                              PnlVect *y, double hx, double ht, PnlVect *V,
                              double Np, double T, double support_espace);
extern void derive_x_operateur_noyau(PnlVect *res, PnlMat *X, PnlVect *T_alea, double s,
                                     PnlVect *y, double hx, double ht, PnlVect *V,
                                     double Np, double T, double support_espace);
extern double derive_t_operateur_noyau( PnlMat *X, PnlVect *T_alea, double s,
                                        PnlVect *y, double hx, double ht, PnlVect *V,
                                        double Np, double T, double support_espace);
extern void derive_xx_operateur_noyau(PnlMat *res, PnlMat *X, PnlVect *T_alea, double s,
                                      PnlVect *y, double hx, double ht, PnlVect *V,
                                      double Np, double T, double support_espace);
#endif /* _NOYAU_BSDE_H */
