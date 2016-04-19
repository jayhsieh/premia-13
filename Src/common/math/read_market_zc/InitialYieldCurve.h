
#ifndef INITIALYILEDCURVE_H_INCLUDED
#define INITIALYILEDCURVE_H_INCLUDED

#include "pnl/pnl_vector.h"

#define INC 1.0e-5 
// Structure where the initial yield curve is saved.
typedef struct ZCMarketData
{
    char *filename; // Name of the file containing P(t,T) when the curve is not flat.
    int FlatOrMarket; // FlatOrMarket=0 if the initial yield curve is flat
    double Rate; // If FlatOrMarket=0, "Rate" is the constant yield of the curve.

    PnlVect* tm; // Vector of the dates
    PnlVect* Pm; // Vector of ZC price for every date tm[i]

    int Nvalue;  // Number of values read in the file.

}ZCMarketData;

/* InitYieldCurve_flag: Flag to decide to read or not ZC bond datas in "initialyields.dat" */
void SetInitYieldCurve(int InitYieldCurve_flag, double R_flat, ZCMarketData* ZCMarket);

// Read the ZC price from the file "initialyield.dat" and put it in the structure "ZCMarket".
void ReadMarketData(ZCMarketData* ZCMarket);

// Compute the ZC price P(0,T) by interpolating the initial yield curve contained in ZCMarket.
double BondPrice(double T, ZCMarketData* ZCMarket);

// Compute f(0, T) the forward rate, known at 0, maturing at T.
double ForwardRate(double T, ZCMarketData* ZCMarket);

// Delete the structure ZCMarket
int DeleteZCMarketData(ZCMarketData* ZCMarket);

// Computes the ATM swaption strike.
double ATMSwaptionStrike(double T_start, double T_end, double period, ZCMarketData* ZCMarket);

#endif /* INITIALYILEDCURVE_H_INCLUDED */

