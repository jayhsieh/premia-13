#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_vector.h"
#include "nelsen_siegel.h"
#include "read_market_data.h"


///******************* ZERO COUPON CURVE *********************///
// Read the ZC price from the file "initialyield.dat" and put it in the structure "ZCMarket".
void ReadMarketData(ZCMarketData* ZCMarket)
{
    FILE* Entrees;                   /*File variable of the code*/

    int i, etat;
    char ligne[100];
    char* pligne;
    double p, tt;

    char data_file[]="brigo_mercurio_zc.dat";
    Entrees=fopen(data_file, "r");

    if(Entrees==NULL)
    {
      printf("Le FICHIER %s N'A PU ETRE OUVERT. VERIFIER LE CHEMIN\n", data_file); abort();
    }

    i=0; // i represents the number of value read in the file
    pligne=ligne;

    ZCMarket->Pm = pnl_vect_create(100);
    ZCMarket->tm = pnl_vect_create(100);
    ZCMarket->NeslenSiegelParams = pnl_vect_create(4);

    while(1)
    {
        pligne=fgets(ligne, sizeof(ligne), Entrees);
        if(pligne==NULL)
        {
            break;
        }
        else
        {
            sscanf(ligne,"%lf t=%lf", &p, &tt);
            /* La ligne lue dans le fichier doit etre de la forme "0.943290 t=0.5" ou 0.943290 est un double pour le prix de B(0,t=0.5)*/
            LET(ZCMarket->Pm,i) = p;   /*enregistre le prix du zero coupon*/
            LET(ZCMarket->tm,i) = tt;  /*enreristre le temps correspondant*/
            i++;
        }
    }

    etat=fclose(Entrees);

    ZCMarket->Nvalue = i;
    pnl_vect_resize(ZCMarket->Pm, i);
    pnl_vect_resize(ZCMarket->tm, i);

    NelsonSiegel_Fitting(ZCMarket);
}

// Compute the ZC price P(0,T) by interpolating the initial yield curve contained in ZCMarket.
double BondPrice(double T, ZCMarketData* ZCMarket)
{
    double beta0, beta1, beta2, tau;

    beta0 = GET(ZCMarket->NeslenSiegelParams, 0);
    beta1 = GET(ZCMarket->NeslenSiegelParams, 1);
    beta2 = GET(ZCMarket->NeslenSiegelParams, 2);
    tau   = GET(ZCMarket->NeslenSiegelParams, 3);

    return BondPrice_NelsonSiegel(T, beta0, beta1, beta2, tau);
}

double ForwardRate(double T, ZCMarketData* ZCMarket)
{
    double eps = 0.01*T;
    return -(log(BondPrice(T + eps,ZCMarket))-log(BondPrice(T,ZCMarket)))/(eps);
}


int DeleteZCMarketData(ZCMarketData* ZCMarket)
{
    if(ZCMarket->FlatOrMarket!=0)
    {
        pnl_vect_free(&(ZCMarket->tm));
        pnl_vect_free(&(ZCMarket->Pm));
        pnl_vect_free(&(ZCMarket->NeslenSiegelParams));
    }

    return 1;
}


///******************* CAP VOLATILITIES *********************///
// Read the caplet volatilities from file and put it in the structure "MktATMCapVolData".
void ReadCapMarketData(MktATMCapVolData *MktATMCapVol)
{
    FILE* Entrees;                   /*File variable of the code*/

    int i, etat;
    char ligne[100];
    char* pligne;
    double p, t0, tn, periodicity;
    char data_file[]="brigo_mercurio_cap.dat";

    Entrees=fopen(data_file, "r");

    if(Entrees==NULL)
    {
      printf("Le FICHIER %s N'A PU ETRE OUVERT. VERIFIER LE CHEMIN\n", data_file); abort();
    }

    i=0; // i represents the number of value read in the file
    pligne=ligne;

    MktATMCapVol->CapVolatility = pnl_vect_create(100);
    MktATMCapVol->CapFirstResetDate = pnl_vect_create(100);
    MktATMCapVol->CapMaturity = pnl_vect_create(100);
    MktATMCapVol->CapPeriodicity = pnl_vect_create(100);

    while(1)
    {
        pligne=fgets(ligne, sizeof(ligne), Entrees);
        if(pligne==NULL)
        {
            break;
        }
        else
        {
            sscanf(ligne, "%lf T0=%lf Tn=%lf periodicity=%lf", &p, &t0, &tn, &periodicity);
            LET(MktATMCapVol->CapVolatility, i) = p; /*Store the cap volatility*/
            LET(MktATMCapVol->CapFirstResetDate, i) = t0;  /*Store the cap maturity*/
            LET(MktATMCapVol->CapMaturity, i) = tn; /*Store the FirstResetDate*/
            LET(MktATMCapVol->CapPeriodicity, i) = periodicity;  /*Store the periodicity*/
            i++;
        }
    }

    etat=fclose(Entrees);

    MktATMCapVol->NbrData = i;

    pnl_vect_resize(MktATMCapVol->CapVolatility, i);
    pnl_vect_resize(MktATMCapVol->CapFirstResetDate, i);
    pnl_vect_resize(MktATMCapVol->CapMaturity, i);
    pnl_vect_resize(MktATMCapVol->CapPeriodicity, i);
}


// Delete caplets data
int DeleteMktATMCapVolData(MktATMCapVolData* MktATMCapVol)
{
    pnl_vect_free(&(MktATMCapVol->CapVolatility));
    pnl_vect_free(&(MktATMCapVol->CapFirstResetDate));
    pnl_vect_free(&(MktATMCapVol->CapMaturity));
    pnl_vect_free(&(MktATMCapVol->CapPeriodicity));

    return 1;
}

///******************* SWAPTION VOLATILITIES *********************///
// Read the Swaption volatilities from file and put it in the structure "MktATMSwaptionVolData".
void ReadSwaptionMarketData(MktATMSwaptionVolData* MktATMSwaptionVol)
{
    FILE* Entrees;                   /*File variable of the code*/

    int i, etat;
    char ligne[100];
    char* pligne;
    double p, t0, tn;

    char data_file[]="brigo_mercurio_swaption.dat";
    Entrees=fopen(data_file, "r");

    if(Entrees==NULL)
    {
      printf("Le FICHIER %s N'A PU ETRE OUVERT. VERIFIER LE CHEMIN\n", data_file); abort();
    }

    i=0; // i represents the number of value read in the file
    pligne=ligne;

    MktATMSwaptionVol->SwaptionMaturity = pnl_vect_create(100);
    MktATMSwaptionVol->SwaptionTenor = pnl_vect_create(100);
    MktATMSwaptionVol->SwaptionVolatility = pnl_vect_create(100);

    pligne=fgets(ligne, sizeof(ligne), Entrees);
    sscanf(ligne, "periodicity=%lf", &t0);
    MktATMSwaptionVol->Periodicity = t0;

    while(1)
    {
        pligne=fgets(ligne, sizeof(ligne), Entrees);
        if(pligne==NULL)
        {
            break;
        }
        else
        {
            sscanf(ligne, "%lf T0=%lf Tn=%lf", &p, &t0, &tn);
            LET(MktATMSwaptionVol->SwaptionVolatility,i) = p; /*Store the Swaption volatility*/
            LET(MktATMSwaptionVol->SwaptionMaturity,i) = t0; /*Store the Swaption maturity*/
            LET(MktATMSwaptionVol->SwaptionTenor,i) = t0+tn;  /*Store the Swaption maturity*/
            //printf("ligne = %s\n", ligne);getchar();
            //printf(" t0 = %f, tn = %f, vol=%f \n", t0, tn, p);getchar();
            i++;
        }
    }

    etat=fclose(Entrees);

    MktATMSwaptionVol->NbrData = i;

    pnl_vect_resize(MktATMSwaptionVol->SwaptionVolatility, i);
    pnl_vect_resize(MktATMSwaptionVol->SwaptionMaturity, i);
    pnl_vect_resize(MktATMSwaptionVol->SwaptionTenor, i);
}

// Delete caplets data
int DeleteMktATMSwaptionVolData(MktATMSwaptionVolData* MktATMSwaptionVol)
{
    pnl_vect_free(&(MktATMSwaptionVol->SwaptionVolatility));
    pnl_vect_free(&(MktATMSwaptionVol->SwaptionMaturity));
    pnl_vect_free(&(MktATMSwaptionVol->SwaptionTenor));

    return 1;
}



/*double BondPrice(double T, ZCMarketData* ZCMarket)
{
    double P0T;
    int i=0;

    if(T>0)
    {
        if(ZCMarket->FlatOrMarket==0) // If there is no curve to read. ie : the initial yield curve is flat.
        {
            P0T = exp(-ZCMarket->Rate * T);
        }
        else
        {
            for(i=0; i<ZCMarket->Nvalue; i++)
            {
                if(T<=GET(ZCMarket->tm,i)) break;
            }

            if(i < ZCMarket->Nvalue)
            {
                P0T = GET(ZCMarket->Pm,i-1)*(GET(ZCMarket->tm,i)-T)/(GET(ZCMarket->tm,i)-GET(ZCMarket->tm,i-1)) + GET(ZCMarket->Pm,i)*(T-GET(ZCMarket->tm,i-1))/(GET(ZCMarket->tm,i)-GET(ZCMarket->tm,i-1));
            }
            else
            {
                P0T=GET(ZCMarket->Pm,i-1)+(T-GET(ZCMarket->tm,i-1))*(GET(ZCMarket->Pm,i-1)-GET(ZCMarket->Pm,i-2))/(GET(ZCMarket->tm,i-1)-GET(ZCMarket->tm,i-2));
            }
        }
    }

    else // P(0,0) = 1
    {
        P0T=1;
    }

    return P0T;
}*/

