
#ifndef TreeLRS1D_H_INCLUDED
#define TreeLRS1D_H_INCLUDED

#include "pnl/pnl_vector.h"
#include "math/read_market_zc/InitialYieldCurve.h"

///*******************TreeLRS1D structure******************///
typedef struct TreeLRS1D
{
  double Tf;            // Final time of the tree, dt=Tf/Ngrid
  int Ngrid;            // Number of time step in the TreeLRS1D

  PnlVect* t;           // Time step grid, from t[0] to T[Ngrid].

  PnlVect *phi;

}TreeLRS1D;

///************* Datas specific to Hull and White **************///
typedef struct ModelLRS1D
{
    double Sigma;
    double Rho;
    double Kappa;
    double Lambda;

}ModelLRS1D;

///******** Fonctions relatives a la construction de l'arbre ********///

int SetTimegridCapLRS1D(TreeLRS1D *Meth, int NtY, double current_date, double T0, double S0, double periodicity);

//Construction of the time grid
int SetTimegridZCbondLRS1D(TreeLRS1D *Meth, int n, double current_date, double T, double S);


// Construction of the time grid
int SetTimegridLRS1D(TreeLRS1D *Meth, int n, double current_date, double T);

void SetTreeLRS1D(TreeLRS1D* Meth, ModelLRS1D* ModelParam, ZCMarketData* ZCMarket);

double r_to_y(ModelLRS1D* ModelParam, double r);

double y_to_r(ModelLRS1D* ModelParam, double y);

/*Compute m, mean of Y=log(r/sigma)*/
double mean(double time,double Y,double Phi, ZCMarketData* ZCMarket, ModelLRS1D* ModelParam);

void probabilities(double date, double y_ij, double phi_ij, double lambda, double sqrt_delta_t, ModelLRS1D* ModelParam, ZCMarketData* ZCMarket, PnlVect* proba_from_ij);

int indice(int i, int h);
double phi_value(TreeLRS1D *Meth, int i, int h, int j); // i>1 , j=0,1,2
double Interpolation(TreeLRS1D *Meth, int i, int h, PnlVect* OptionPriceVect2, double phi_star);
double MeanPrice(TreeLRS1D *Meth, int i, int h, PnlVect* OptionPriceVect2);
int number_phi_in_box(int i, int h);
int index_tree(int i, int h, int j);

int indiceTimeLRS1D(TreeLRS1D *Meth, double s); // To locate the date s inf the tree. t[indiceTimeLRS1D(s)-1]< s <= t[indiceTimeLRS1D(s)]

int DeleteTimegridLRS1D(struct TreeLRS1D *Meth);

int DeleteTreeLRS1D(struct TreeLRS1D* Meth);

#endif // HW2DTREE_H_INCLUDED
