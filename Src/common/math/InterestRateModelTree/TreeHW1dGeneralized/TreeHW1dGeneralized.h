
#ifndef TreeHW1dGeneralized_H_INCLUDED
#define TreeHW1dGeneralized_H_INCLUDED
#include "math/read_market_zc/InitialYieldCurve.h"



///*******************TreeHW1dG structure******************///
typedef struct TreeHW1dG
{
  double Tf;            // Final time of the tree, dt=Tf/Ngrid
  int Ngrid;            // Number of time step in the TreeHW1dG

  PnlVect* t;           // Time step grid, from t[0] to T[Ngrid].
  PnlVectInt* Jminimum; // Jminimum[i] : Minimal index at time i
  PnlVectInt* Jmaximum; // Jmaximum[i] : Maximal index at time i
  PnlVect* alpha;       // Translation from x to r. ( r_t = x_t + alpha_t)
}TreeHW1dG;

///************* Datas specific to Hull and White **************///
typedef struct ModelHW1dG
{
    double MeanReversion;                     /*Speed revertion of the Hullwhite model.*/

    PnlVect* TimeGrid;
    PnlVect* ShortRateVolGrid; /*Volatility of the Hullwhite model.*/

}ModelHW1dG;

double Current_VolatilityHW1dG(ModelHW1dG* HW1dG, double t);

void SetTreeHW1dG(TreeHW1dG* Meth, ModelHW1dG* ModelParam, ZCMarketData* ZCMarket);

int SetTimeGridHW1dG(TreeHW1dG *Meth, int NbrTimeStep, double T1, double T2);

int SetTimeGrid_TenorHW1dG(TreeHW1dG *Meth, int NtY, double T0, double S0, double periodicity);

void BackwardIterationHW1dG(TreeHW1dG* Meth, ModelHW1dG* ModelParam, PnlVect* OptionPriceVect1, PnlVect* OptionPriceVect2, int index_last, int index_first);

double SpaceStepHW1dG(double delta_t, double sigma); // Renvoie Delta_x(i)

double ProbaUpHW1dG(int j, int k, double delta_t2,double beta_x, double mean_reversion); // beta_x = deltax1/deltax2

double ProbaMiddleHW1dG(int j, int k, double delta_t2,double beta_x, double mean_reversion);

double ProbaDownHW1dG(int j, int k, double delta_t2, double beta_x, double mean_reversion);

int IndexTimeHW1dG(TreeHW1dG *Meth, double s); // To locate the date s inf the tree.

int DeleteTimegridHW1dG(struct TreeHW1dG *Meth);

int DeleteTreeHW1dG(struct TreeHW1dG* Meth);

int DeletModelHW1dG(struct ModelHW1dG* HW1dG);


#endif // TreeHW1dGeneralized_H_INCLUDED
