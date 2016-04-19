#ifndef READ_MARKET_DATA_H_INCLUDED
#define READ_MARKET_DATA_H_INCLUDED

#include "pnl/pnl_vector.h"

// Structure where cap volatilities is saved.
typedef struct MktATMCapVolData
{
    PnlVect* CapPeriodicity;
    PnlVect* CapFirstResetDate; // Vector of the CapFirstResetDate
    PnlVect* CapMaturity; // Vector of the maturities
    PnlVect* CapVolatility; // Vector of Cap Volatilities for every maturity CapMaturity[i]

    int NbrData;  // Number of values read in the file.

}MktATMCapVolData;

// Structure where swaption volatilities is saved.
typedef struct MktATMSwaptionVolData
{
    double Periodicity;
    PnlVect* SwaptionMaturity;    // Vector of the maturities
    PnlVect* SwaptionTenor;    // Vector of the maturities
    PnlVect* SwaptionVolatility;  // Vector of Cap Volatilities for every maturity CapMaturity[i]

    int NbrData;  // Number of values read in the file.

}MktATMSwaptionVolData;


// Structure where the initial yield curve is saved.
typedef struct ZCMarketData
{
    int FlatOrMarket; // FlatOrMarket=0 if the initial yield curve is flat
    double Rate; // If FlatOrMarket=0, "Rate" is the constant yield of the curve.

    PnlVect* tm; // Vector of the dates
    PnlVect* Pm; // Vector of ZC price for every date tm[i]
    PnlVect* NeslenSiegelParams;

    int Nvalue;  // Number of values read in the file.

}ZCMarketData;

// Structure to be use when calibrating on swaptions
typedef struct {
  ZCMarketData *zc_market;
  MktATMSwaptionVolData *swaption_vol_market;
} MktSwpZCData;

// Structure to be use when calibrating on caps
typedef struct {
  ZCMarketData *zc_market;
  MktATMCapVolData *cap_vol_market;
} MktCapZCData;

///********* Read zero coupon curve from file ***********///
// Read ZC prices from file and store it in a structure "ZCMarketData".
void ReadMarketData(ZCMarketData* ZCMarket);
// Compute the ZC price P(0,T) by interpolating the initial yield curve contained in ZCMarket.
double BondPrice(double T, ZCMarketData* ZCMarket);
int DeleteZCMarketData(ZCMarketData* ZCMarket);
double ForwardRate(double T, ZCMarketData* ZCMarket);

///********* Read cap volatilities from file ***********///
// Read cap volatilities from file and store it in a structure "MktATMCapVolData".
void ReadCapMarketData(MktATMCapVolData* MktATMCapVol);
int DeleteMktATMCapVolData(MktATMCapVolData* MktATMCapVol);

///********* Read swaption volatilities from file ***********///
// Read swaption volatilities from file and store it in a structure "MktATMSwaptionVolData".
void ReadSwaptionMarketData(MktATMSwaptionVolData* MktATMSwaptionVol);
int DeleteMktATMSwaptionVolData(MktATMSwaptionVolData* MktATMSwaptionVol);

#endif // HW1DGCALIBRATION_H_INCLUDED
