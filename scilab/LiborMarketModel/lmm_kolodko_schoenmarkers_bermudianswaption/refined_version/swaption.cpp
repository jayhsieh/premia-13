/*---------------------------------------*/
/*   BERMUDAN SWAPTION PRICER            */
/*    Kolodko & Schoenmakers Algorithm   */
/*    Andersen Algorithm                 */
/*    Tight upper estimator              */
/*                                       */
/*---------------------------------------*/
/*  Julien Bourgouint, Premia July 2006  */
/*  julien_bourgouint@yahoo.fr           */ 
/*---------------------------------------*/

#include <iostream>
#include <fstream>
#include <valarray>
using namespace std;
#include "swaption.h"
#include "proba.h"

LMM::LMM(int _first, int _k, double _L_0,  int _d, double phi, double _a, double _b, double _c, double _ginf, double _delta, int _p, double _theta):first(_first), kinit(_k),d(_d),M(_k,phi,_d), a(_a), b(_b), c(_c), ginf(_ginf), delta(_delta), p(_p), theta(_theta), dt(delta/p), L_0(_L_0), sqrt_dt(sqrt(dt))
/* constructor for LIBOR rates */
{
  L.resize(_k);
  M.reduce();
  /* gamma coefficients at different times*/
  double s;
  Gamma.resize((_k-1)*p);
  for (int j = 0; j < (_k-1)*p; j++)
    {
      s = j*dt;
      Gamma[j] = c*(ginf + (1 - ginf+ a*s)*exp(-b*s));
    }
  /* initialization of actualization coefficients */
  actual.resize(_k+1);
  actual[0] = 1.;  
  for (int i = 1; i <= _k; i++) 
    {
      actual[i] = actual[i-1]/(1+delta*L_0);
    }
  B.resize(_k+1);  B[0] = 1.;
  exp_b.resize(2*_k);  
  exp_b[0] = 1.;
  double eb = exp(delta*b);
  for (int i = 1; i < 2*_k; i++)
    { 
      exp_b[i] = exp_b[i - 1]*eb;
    }
};

Y1::Y1(int first, int k, double L_0, int d, double phi, double a, double b, double c, double ginf, double delta, int p, double theta, int _N):LMM(first, k, L_0, d, phi, a, b, c, ginf, delta, p, theta), N(_N)
{
};

Y2::Y2 (int first, int k, double L_0, int d, double phi, double a, double b, double c, double ginf, double delta, int p, double theta, int _N1, int _N2):LMM(first, k, L_0, d, phi, a, b, c, ginf, delta, p, theta), N1(_N1), N2(_N2)
{
  P.resize(k);	
  Lmem.resize(k);
};

Andersen::Andersen(int _strat, int first, int k, double L_0, int d, double phi, double a, double b, double c, double ginf, double delta, int p, double theta, int _div, int _NH, int _NS):strat(_strat), LMM(first, k, L_0, d, phi, a, b, c, ginf, delta, p, theta), div(_div), NH(_NH), NS(_NS)
{  
  H.resize(k); H = .1; H[k-1] = 0.;// H is initialized with a high value (>> payoff)
};


Andersenhat::Andersenhat(int _N1, int _N2, int strat, int first, int k, double L_0, int d, double phi, double a, double b, double c, double ginf, double delta, int p, double theta, int div, int NH):N1(_N1),N2(_N2), Andersen(strat, first, k, L_0, d, phi, a, b, c, ginf, delta, p, theta, div, NH, 0)
{
  P.resize(k);	
  Lmem.resize(k);
};

Y1up::Y1up(int first, int k, double L_0, int d, double phi, double a, double b, double c, double ginf, double delta, int p, double theta, int _NK, int _N1up, int _N2up):LMM(first, k, L_0, d, phi, a, b, c, ginf, delta, p, theta), NK(_NK), N1up(_N1up), N2up(_N2up)
{
  Lmem1.resize(k);
  Lmem2.resize(k);
};

YAup::YAup(int strat, int first, int k, double L_0, int d, double phi, double a, double b, double c, double ginf, double delta, int p, double theta, int div, int NH, int _NK, int _N1up, int _N2up):Andersen(strat, first, k, L_0, d, phi, a, b, c, ginf, delta, p, theta, div, NH, 0), NK(_NK), N1up(_N1up), N2up(_N2up)
{
  Lmem1.resize(k);
  Lmem2.resize(k);
};

void LMM::Scheme()
/* log Euler Scheme to compute LIBOR's rates between T_i and T_i+1 */
{ 
  k--;//remaining exercice dates
  int l,i,j;   
  /* copying of the still alive LIBOR rates in the L first boxes */   
  for (i = 0; i < k; i++)
	{  
	  L[i] = L[i+1];
	}
  /* rates dynamic */ 
  double brown[d]; 
  double det, dw;  
  for (l = 0; l < p; l++) //time path
    {   
      for (i = 0; i < d; i++)
	{ 
	  brown[i] = sqrt_dt*G();
	}
      for (i = k-1; i >= 0; i--)
	{ 
	  det = 0.;
	  for (j = 0; j <= i; j++)
	    {
	      det += L[j]*Gamma[(i+1)*p-l]*M.ro_d[i][j]*Gamma[(j+1)*p - l]/(1+delta*L[j]);
	    }
	  det *= delta;
	  dw = 0.;
	  for (j = 0; j < d; j++)
	    { 
	      dw += brown[j]*M.les_ei[i][j];
	    }
	  L[i] *= 1 + det*dt + Gamma[(i+1)*p-l]*dw;
	}
    } 
  /* payoff computation */
  payoff = 0.;
  if (k <= kinit -  first)
    {
      for (j = 0; j < k; j++)
	{
	  B[j+1] = B[j]/(1+delta*L[j]);
	  payoff += B[j+1]*(L[j]-theta);	
	}
      payoff = max(payoff, 0.)*delta*actual[kinit - k];//actualized pay off
    }
}

double LMM::CF()
/* 
Using the Black approximation formula for european swaption prices in the LIBOR market model, 
we return the highest european price for maturities ranging from current date to 
terminal time.  
*/
{
  int pp; //named pp to avoid confusion with p (precision of the Euler Scheme)
  int i, j;
  double Bpq, Spq, sigmapq_carre, integrale, sigmapq, dplus, dmoins, result;
  double maxi = 0.;
  for (pp = 1; pp < k; pp++)
    { 
      Bpq = 0.;
      for (i = pp; i < k; i++) Bpq += delta*B[i+1]; 
      Spq = (B[pp] - B[k])/Bpq; 
      sigmapq_carre =  0.;   
      for (i = pp; i < k; i++) 
	{
	  for (j = pp; j < k; j++) 
	    {
	      /* computed formula */
	      double Ti = i*delta;
	      double Tj = j*delta;
	      double Tp = pp*delta;
	      double Tm = (Ti+Tj)/2.;
	      double DT = abs(Ti - Tj)/2.; 
	      integrale = ginf*ginf*Tp +  ginf/(exp_b[i]*b)*( (1-ginf + a/b)*(exp_b[pp] - 1) - a*((Tp-Ti)*exp_b[pp] - Ti) ) + ginf/(exp_b[j]*b)*( (1-ginf + a/b)*(exp_b[pp] - 1) - a*((Tp-Tj)*exp_b[pp] - Tj) ) + 1/(exp_b[i+j])* ( (1-ginf)*(1-ginf)*(exp_b[2*pp] - 1)/(2*b) + a*(1-ginf)/b*( - ((Tp- Tm)* exp_b[2*pp] + Tm) + (exp_b[2*pp] - 1)/(2*b)) + a*a*DT*DT/(2*b)*(exp_b[2*pp] - 1) +  a*a/(8*b*b*b)*(exp_b[2*pp]*(4*b*b*(Tp - Tm)*(Tp-Tm) - 2*2*b*(Tp-Tm) +2 ) - 4*b*b*Tm*Tm + 4*b*Tm + 2));
	      integrale *=c*c;
	      //trapezoidale rule (less efficient altenative for the computed formula)
	      /*integrale = 0.;
	      for (int time = 0; time < pp*p; time++)
		{
		integrale += (Gamma[j*p - time]*Gamma[i*p - time]+Gamma[j*p - time-1]*Gamma[i*p - time-1])/2.;
		}
              integrale *= delta/p;*/
	      integrale *= M.ro_d[j-1][i-1]; 
	      sigmapq_carre += B[i]*B[j]*L[i-1]*L[j-1]/((B[pp] - B[k])*(B[pp] - B[k]))*integrale;
	    }
	}
      sigmapq = sqrt(sigmapq_carre)*delta; 
      dplus = log(Spq/theta)/sigmapq + sigmapq/2;
      dmoins = dplus - sigmapq;
      result = Bpq*(Spq*Cumul_Normal(dplus) - theta*Cumul_Normal(dmoins));
      maxi = max(maxi, result);
    }
  return (maxi*actual[kinit-k]);
}

double LMM::next_euro()
/* 
Using the Black approximation formula for european swaption prices in the LIBOR market model, 
we return the european price with maturity equal to the current date + delta. 
*/
{
  int pp = 1; //so named to avoid confusion with p (precision of the Euler Scheme)
  int i, j;
  double Bpq, Spq, sigmapq_carre, integrale, sigmapq, dplus, dmoins, result;
  double maxi = 0.;
  Bpq = 0.;
  for (i = pp; i < k; i++) Bpq += delta*B[i+1]; 
  Spq = (B[pp] - B[k])/Bpq; 
  sigmapq_carre =  0.;   
  for (i = pp; i < k; i++) 
    {
      for (j = pp; j < k; j++) 
	{
	  /* computed formula */
	  double Ti = i*delta;
	  double Tj = j*delta;
	  double Tp = pp*delta;
	  double Tm = (Ti+Tj)/2.;
	  double DT = abs(Ti - Tj)/2.; 
	  integrale = ginf*ginf*Tp +  ginf/(exp_b[i]*b)*( (1-ginf + a/b)*(exp_b[pp] - 1) - a*((Tp-Ti)*exp_b[pp] - Ti) ) + ginf/(exp_b[j]*b)*( (1-ginf + a/b)*(exp_b[pp] - 1) - a*((Tp-Tj)*exp_b[pp] - Tj) ) + 1/(exp_b[i+j])* ( (1-ginf)*(1-ginf)*(exp_b[2*pp] - 1)/(2*b) + a*(1-ginf)/b*( - ((Tp- Tm)* exp_b[2*pp] + Tm) + (exp_b[2*pp] - 1)/(2*b)) + a*a*DT*DT/(2*b)*(exp_b[2*pp] - 1) +  a*a/(8*b*b*b)*(exp_b[2*pp]*(4*b*b*(Tp - Tm)*(Tp-Tm) - 2*2*b*(Tp-Tm) +2 ) - 4*b*b*Tm*Tm + 4*b*Tm + 2));
	  integrale *=c*c;
	  //trapezoidale rule (less efficient altenative for the computed formula)
	  /*integrale = 0.;
	    for (int time = 0; time < pp*p; time++)
	    {
	    integrale += (Gamma[j*p - time]*Gamma[i*p - time]+Gamma[j*p - time-1]*Gamma[i*p - time-1])/2.;
	    }
	    integrale *= delta/p;*/
	  integrale *= M.ro_d[j-1][i-1]; 
	  sigmapq_carre += B[i]*B[j]*L[i-1]*L[j-1]/((B[pp] - B[k])*(B[pp] - B[k]))*integrale;
	}
    }
  sigmapq = sqrt(sigmapq_carre)*delta; 
  dplus = log(Spq/theta)/sigmapq + sigmapq/2;
  dmoins = dplus - sigmapq;
  result = Bpq*(Spq*Cumul_Normal(dplus) - theta*Cumul_Normal(dmoins));
  return (result*actual[kinit-k]);
}

void LMM::InitialCond()
{
  k = kinit;
  L = L_0;
  /* payoff computing */ 
  payoff = 0.;
  for (int j = 0; j < k; j++)
    {
      B[j+1] = B[j]/(1+delta*L[j]);
      payoff += B[j+1]*(L[j]-theta);	
    }
  payoff = max(payoff, 0.)*delta*actual[kinit - k];//actualised pay off
}

void LMM::print_CF()
/* 
Using the Black approximation formula for european swaption prices in the LIBOR market model, 
we compute and print out the highest european price for maturities ranging from current date to 
terminal time and create a scilab file CF.sci with the table of european prices.
*/
{
  int pp, i, j, temps;
  double Bpq, Spq, sigmapq_carre, integrale, sigmapq, dplus, dmoins, result;
  for (int j = 0; j < k; j++)
    {
      B[j+1] = B[j]/(1+delta*L[j]);
    }
  valarray<double> tab;
  tab.resize(k - 1); tab = 0.;
  for (pp = 1; pp < k; pp++)
    { 
      Bpq = 0.;
      for (i = pp; i < k; i++) Bpq += delta*B[i+1]; 
      Spq = (B[pp] - B[k])/Bpq; 
      sigmapq_carre =  0.;   
      for (i = pp; i < k; i++) 
	{
	  for (j = pp; j < k; j++) 
	    {
	      /* computed formula */ 

	      double Ti = i*delta;
	      double Tj = j*delta;
	      double Tp = pp*delta;
	      double Tm = (Ti+Tj)/2.;
	      double DT = abs(Ti - Tj)/2.; 
	      integrale = ginf*ginf*Tp +  ginf/(exp_b[i]*b)*( (1-ginf + a/b)*(exp_b[pp] - 1) - a*((Tp-Ti)*exp_b[pp] - Ti) ) + ginf/(exp_b[j]*b)*( (1-ginf + a/b)*(exp_b[pp] - 1) - a*((Tp-Tj)*exp_b[pp] - Tj) ) + 1/(exp_b[i+j])* ( (1-ginf)*(1-ginf)*(exp_b[2*pp] - 1)/(2*b) + a*(1-ginf)/b*( - ((Tp- Tm)* exp_b[2*pp] + Tm) + (exp_b[2*pp] - 1)/(2*b)) + a*a*DT*DT/(2*b)*(exp_b[2*pp] - 1) +  a*a/(8*b*b*b)*(exp_b[2*pp]*(4*b*b*(Tp - Tm)*(Tp-Tm) - 2*2*b*(Tp-Tm) +2 ) - 4*b*b*Tm*Tm + 4*b*Tm + 2));
	      integrale *=c*c;
	      //trapezoidale rule
	      /*integrale = 0.;
	      for (temps = 0; temps < pp*p; temps++)
		{
		 integrale += (Gamma[j*p - temps]*Gamma[i*p - temps]+Gamma[j*p - temps-1]*Gamma[i*p - temps-1])/2.;
		}
	      integrale *= delta/p*/
	      integrale *= M.ro_d[j-1][i-1];
	      sigmapq_carre += B[i]*B[j]*L[i-1]*L[j-1]/((B[pp] - B[k])*(B[pp] - B[k]))*integrale;
	    }
	}
      sigmapq = sqrt(sigmapq_carre)*delta;
      dplus = log(Spq/theta)/sigmapq + sigmapq/2;
      dmoins = dplus - sigmapq;
      tab[pp - 1] = Bpq*(Spq*Cumul_Normal(dplus) - theta*Cumul_Normal(dmoins));
    }
  tab *= 10000; //in base points
  ofstream monfichier("CF.sci");   
  monfichier << "function A = CF()" << endl << "A = ["; 
  for (pp = 0; pp < k - 1; pp++)
    { 
      monfichier << tab[pp];
      if (pp != k - 2) monfichier << ", ";
    }
  monfichier << "];" << endl << " return A" << endl << "endfunction";
  cout << "lower bound : best european (Black formula) : " << tab.max() << endl;
}

void LMM::print_MC(int N)
/* 
Using Monte Carlo simulations, we compute and print out the highest european price for 
maturities ranging from current date to terminal time and create a scilab file MC.sci 
with the table of european prices. 
*/
{
  valarray<double> tab;
  tab.resize(k); tab = 0.;
  int kref = k;
  for (int j = 0; j < N; j++)
    {
      InitialCond();
      for(int i = 0; i < kref; i++)
	{ 
	  Scheme();
	  tab[i] += payoff;
	}
    }
  tab /= N;
  tab *= 10000; //in base points
  ofstream monfichier("MC.sci");   
  monfichier << "function A = MC()" << endl << "A = ["; 
  for(int i = 0; i < kref; i++)
    { 
      monfichier << tab[i];
      if (i != kref-1) monfichier << ", ";
    }
  monfichier << "];" << endl << " return A" << endl << "endfunction";
  cout << "lower bound : best european (Monte Carlo) : " << tab.max() << endl;
}

double LMM::trajup()
{
/* computes one rough upper bound trajectory */
  InitialCond();
  double maxi = 0.;
  int kref = k;
  for(int i = 0; i < kref; i++)
    { 
      Scheme();
      maxi = max(payoff, maxi);
    }
  return maxi;
}

void LMM::estimateup(int N)
/* rough upper bound estimation (provide arithmetic mean, standard deviation and a 95% confidence interval) */
{
  double traj_up;
  double mean = 0.;
  double var = 0.;
  for (int i = 0; i < N; i++) 
    {
      traj_up = trajup();
      mean += traj_up;
      var += traj_up*traj_up;
    }  
  mean /= N;
  var /= N;
  double SD = sqrt(var-mean*mean);
  double inter = 1.96*SD/sqrt(N);
  inter *= 10000; mean *= 10000;//result is in base points
  cout << "rough upper bound :" << mean << "  (95% confidence interval: [" << (mean - inter) << "," << (mean + inter) << "]" << " )" << endl;
}

double Y1::short_traj1()
/* simulates Z_tau1 */
{
  bool b = true; //false = exercice
  for(int i = 0; (i < kinit - 1 ) && b; i++)
    {    
      Scheme();  
      if (payoff != 0.) //no need to compare  the payoff in this case (frequent in OTM case)
	{
	  b = (payoff <= CF());
	}
    }     
  return payoff;
}

void Y1::estimate()
/* provides arithmetic mean, standard deviation and a 95% confidence interval for the price estimator */
{
  double traj;
  double mean = 0.;
  double var = 0.;
  for(int i = 0; i < N; i++)
    {
      InitialCond(); 
      traj = short_traj1();
      mean += traj;    
      var += traj*traj;
    }
  mean /= N; var /= N;
  var -= mean*mean;
  double SD = sqrt(var);
  double inter = 1.96*SD/sqrt(N); 
  inter *= 10000; mean *= 10000;//result is in base points
  cout << "price tight lower estimator (after 1 iteration) : " << mean << "  (95 % confidence interval : [" << (mean - inter) << "," << (mean + inter) << "]" << " )" << endl;
}

void Y2::traj1()
/* builds a table with Z_tau1 values */
{
  int i; 
  double best_euro;
  /* saves the current rate data */
  int kref = k;
  for(int i = 0; i < k; i++) Lmem[i] = L[i];
  bool A[k-1]; // A[i] false means exercice at time i. table A's length is  k-1 because we exerce in k in all cases.
  for(i = 0; i < kref; i++)
    {    
      Scheme();
      P[i] = payoff; 
      if (i != kref - 1)
	{   
	  best_euro = CF(); 
	  A[i] = (P[i] <= best_euro);
	}
    }     
  for(i = kref - 2; i >= 0; i--)
    { 
      if (A[i]) 
	{ 
	  P[i] = P[i+1];
	} 
    }
  /*  loads the current rate data */
  k = kref;
  for(int i = 0; i < k; i++) L[i] = Lmem[i];
}

double Y2::best_y1()
/* builds a table with values of Y1 and returns their maximum */
{
  valarray<double> tab;
  tab.resize(k); tab = 0.;
  for (int j = 0; j < N2; j++)
    {
      traj1();
      for (int i = 0; i < k; i++)
	{
	  tab[i] += P[i];
	}
    }
  tab /= N2;
  return tab.max();
}

double Y2::short_trajZ()
/* simulates Z_tau2 - Z_tau1 */
{ 
  double payoff1 = 0.;// payoff first iteration
  double po;// payoff
  bool b = true;
  for(int i = 0; (i <= kinit - 2) && b; i++)
    {  
      Scheme();   
      po = payoff;
      if (po != 0.) //no need to compare  the payoff in this case (frequent in OTM case)
	{
	  if (po > CF()) {if (payoff1 == 0.) payoff1 = po;} 
	  b = (po < best_y1());
	}
    }
  return (po - payoff1);
}

void Y2::estimate()
/* provides the arithmetic mean, standard deviation and a 95% confidence interval for the price estimator */
{
  double traj;
  double mean = 0.;
  double var = 0.; 
  for(int i = 0; i < N1; i++)
    {  
      InitialCond(); 
      traj = short_trajZ();
      mean += traj;    
      var += traj*traj;
    }
  mean /= N1; var /= N1;
  var -= mean*mean;
  double SD = sqrt(var);
  double inter = 1.96*SD/sqrt(N1); 
  inter *= 10000; mean *= 10000;//result is in base points
  cout << "(this  following information is more relevant in the refined case)" << endl << "price lower estimator (2 iterations) - price lower estimator (1 iteration) : " <<  mean << "  (95 % confidence interval : [" << (mean - inter) << "," << (mean + inter) << "]" << " )" << endl;
}

double Andersen::Evaluate_H()
/* estimates the price with the current boundary H */
{
  double mean = 0.;
  for (int i = 0; i < NH; i++)
    {
      int j = 0;
      if (strat == 1) {while((j < kinit - 1) && (po[i][j] <= H[j])) j++;}  
      if (strat == 2 or strat == 4) {while((j < kinit - 1) && ((po[i][j] <= H[j]) or (po[i][j] <= euro[i][j]))) j++;}
      if (strat == 3 or strat == 5) {while((j < kinit - 1) && (po[i][j] <= H[j] + euro[i][j])) j++;}
      mean += po[i][j];
    }
  mean /= NH;
  return mean;
}

void Andersen::tab_payoff(int N)
/* builds a payoff table po, each line i corresponds to a trajectory (column j = Time T j) */
{
  if (strat == 1)
    {
      po.resize(N);
      for(int i = 0; i < N; i++)
	{ 
	  po[i].resize(kinit);
	  InitialCond();
	  for(int j = 0; j < kinit; j ++)
	    {
	      Scheme();
	      po[i][j] = payoff; 
	    }
	}  
    }
  if (strat == 2 or strat == 3)
    {
      po.resize(N);
      euro.resize(N);// contains the highest european price for maturities ranging from current date to terminal time.
      for(int i = 0; i < N; i++)
	{ 
	  po[i].resize(kinit);
	  euro[i].resize(kinit);
	  InitialCond();
	  for(int j = 0; j < kinit; j ++)
	    {
	      Scheme();
	      po[i][j] = payoff; 
	      euro[i][j] = CF();
	    }
	}  
    }  
  if (strat == 4 or strat == 5)
    {
      po.resize(N);
      euro.resize(N); // contains the european price for maturity equal current date + delta
      for(int i = 0; i < N; i++)
	{ 
	  po[i].resize(kinit);
	  euro[i].resize(kinit);
	  InitialCond();
	  for(int j = 0; j < kinit; j ++)
	    {
	      Scheme();
	      po[i][j] = payoff; 
	      euro[i][j] = next_euro();
	    }
	}  
    }
}

void Andersen::Move_H(int step, double h)
/* shifts the boundary H keeping its piecewise linearity */
{
  int begin = max(step*(kinit-first + 1)/div - 1, 0) + first - 1;
  H[begin] = h;
  int end = (step+1)*(kinit-first + 1)/div - 1 + first - 1;
  int d = end - begin;
  int i, j;
  for (i = begin+1; i < end; i++)
    {
      H[i] = (H[end]*(i-begin)+H[begin]*(end-i))/d;
    }
}

/*
this golden search of the maximum is adapted from numerical recipes
*/

#define R 0.61803399 //The golden ratios.
#define C (1.0-R)
#define SHFT2(a,b,c) (a)=(b);(b)=(c)
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)

double Andersen::golden(double ax, double tol, int step)//tol is the tolerance for the precision of the result
{
  /* rough bracketing of the minimum */
  double dx = .01;
  double cx = ax + dx;
  double bx = cx + dx;
  Move_H(step, bx);
  double fbx = Evaluate_H();;
  Move_H(step, cx);
  double  fcx = Evaluate_H();
  while(fbx > fcx)
    {
      SHFT3(ax, bx, cx, cx +dx);
      Move_H(step, cx);
      fcx= Evaluate_H();
    }
  /* golden section */
  double f1,f2,x0,x1,x2,x3, fx2, fx1;  
  x0=ax; 
  x3=cx;
  if (fabs(cx-bx) > fabs(bx-ax)) 
    { 
      x1=bx;
      x2=bx+ C*(cx-bx); 
    } 
  else 
    {
      x2=bx;
      x1=bx-C*(bx-ax);
    }
  Move_H(step, x1);
  f1= Evaluate_H();
  Move_H(step, x2);
  f2= Evaluate_H();
  while (fabs(x1 - x2) > tol)  // this criterion may be changed
    {
      if (f2 > f1) 
	{
	  SHFT3(x0,x1,x2,R*x1+C*x3);
	  Move_H(step, x2);
	  fx2 = Evaluate_H();
	  SHFT2(f1,f2,fx2); 
	} 
      else 
	{ 
	  SHFT3(x3,x2,x1,R*x2+C*x0);
	  Move_H(step, x1);
	  fx1 = Evaluate_H();
	  SHFT2(f2,f1,fx1);
	}
    } 
  if (f1 > f2) 
    { 
      return x1;
    } 
  else 
    {
      return x2;
    }
}

void Andersen::Find_H()
/* find the optimal frontier exercice */
{
  tab_payoff(NH);
  double h;
  double ax = 0.;
  for (int step = div-1; step >= 0; step--) // backward determination of the optimal boundary 
    {  
      h = golden(ax, .000001, step);
      ax = h;
      Move_H(step, h);
    } 
}

void Andersen::estimate()
/* provides the arithmetic mean, standard deviation and a 95% confidence interval for the price estimator and the values of the exercise boundary */
{
  Find_H();
  tab_payoff(NS); // an other MC simulation is required (N is taken bigger for the occasion)
  double mean = 0., var = 0.;
  for (int i = 0; i < NS; i++)
    {
      int j = 0;
      if (strat == 1) {while((j < kinit - 1) && (po[i][j] <= H[j])) j++;}
      if (strat == 2 or strat == 4) {while((j < kinit - 1) && ((po[i][j] <= H[j]) or (po[i][j] <= euro[i][j]))) j++;}
      if (strat == 3 or strat == 5) {while((j < kinit - 1) && (po[i][j] <= H[j] + euro[i][j])) j++;}
      mean += po[i][j];
      var += po[i][j]*po[i][j];
    }
  mean /= NS;
  var /= NS;
  double SD = sqrt(var);
  mean *= 10000; SD *= 10000;//result is in base points
  double inter = SD/sqrt(NS);
  cout << "price (Andersen, strategy : " << strat << " ) : " << mean << "  (95 % confidence interval : [" << mean - inter << " , " << mean + inter << " ] " << endl;  
  ofstream monfichier("B.sci");   
  monfichier << "function A = B()" << endl << "A = ["; 
  for (int i = 0; i < kinit; i++)
    { 
      monfichier << 10000*H[i];
      if (i != kinit - 1) monfichier << ", ";
    }
  monfichier << "];" << endl << " return A" << endl << "endfunction";
  cout << endl << "# # # # # # " << endl << "You can get a graph of the exercice boundary by running scilab and writting" << endl << "getf(~/B.sci)" << endl << "plot2d(B())" << endl << "# # # # # # " << endl << endl;
}

void Andersenhat::trajA()
/* builds a table of the payoffs with Andersen exercise strategy */
{
  int i; 
  double best_A;
  /* saves the information */
  int kref = k;
  for(int i = 0; i < k; i++) Lmem[i] = L[i];
  bool A[k-1]; // A[i] false = exercice at time i.
  for(i = 0; i < kref; i++)
    {    
      Scheme();
      P[i] = payoff; 
      if (i != kref - 1)
	{    
	  A[i] = (P[i] <= H[i]);
	}
    }     
  for(i = kref - 2; i >= 0; i--)
    { 
      if (A[i]) 
	{ 
	  P[i] = P[i+1];
	} 
    }
  /*  loads */
  k = kref;
  for(int i = 0; i < k; i++) L[i] = Lmem[i];
}

double Andersenhat::best_yA()
/* returns max(Y,A) */
{
  valarray<double> tab;
  tab.resize(k); tab = 0.;
  for (int j = 0; j < N2; j++)
    {
      trajA();
      for (int i = 0; i < k; i++)
	{
	  tab[i] += P[i];
	}
    }
  tab /= N2;
  return tab.max();
}

double Andersenhat::short_trajZ()
{ 
  double payoffA = 0.;// payoff with Andersen algorithm
  double payoffAhat = 0.;// payoff with Andersen algorithm and 1 iteration of K&S algorithm
  double payo;// payoff
  for(int i = 0; (i <= kinit - 2); i++)
    {  
      Scheme();   
      payo = payoff;
      if (payo != 0.)//no need to compare  the payoff in this case (frequent in OTM case)
	{
	  if (payoffA == 0.) {if (payo > H[i]) payoffA = payo;} 
	  if (payoffAhat == 0.) {if (payo > best_yA()) payoffAhat = payo;}
	}
    }
  return (payoffAhat - payoffA);
}

void Andersenhat::estimate2()
/* provides the arithmetic mean, standard deviation and a 95% confidence interval for the price estimator and the values of the exercice boundary */
{
  Find_H();
  cout << endl;
  double traj;
  double mean = 0.; 
  double var = 0.;
  for(int i = 0; i < N1; i++)
    {  
      InitialCond(); 
      traj = short_trajZ();
      mean += traj;    
      var += traj*traj;
    }
  mean /= N1; var /= N1;
  var -= mean*mean;
  double SD = sqrt(var);
  double inter = 1.96*SD/sqrt(N1); 
  inter *= 10000; mean *= 10000;
  cout << "(this  following information is relevant only in the refined case)" << endl << "price (Andersen, strategy 1 and 1 iteration of K & S) - price (Andersen): " <<  mean << "  (95 % confidence interval : [" << (mean - inter) << "," << (mean + inter) << "]" << " )" << endl;
}

double Y1up::short_traj1()
/* simulates Z_tau1 */
{ 
  /* saves */
  int kref = k;
  for(int i = 0; i < k; i++) Lmem2[i] = L[i];
  bool b = true; // false : exercice
  for(int i = 0; (i < kref -1 ) && b; i++)
    {     
      Scheme(); 
      if (payoff != 0.) // no need to compare  the payoff in this case (frequent in OTM case)
	{
	  b = (payoff <= CF());
	}
    }     
  /* load */
  k = kref;
  for(int i = 0; i < k; i++) L[i] = Lmem2[i];
  return payoff;
}

double Y1up::e_v_short_traj1()
/* calculates an estimation of the price with K&S's algorithm (1 iteration) */
{
 if (payoff > CF()) 
    {
      return payoff;
    }
 else
   {
     double mean = 0.;
     for (int i = 0; i < N2up; i++)
       {
	 mean += short_traj1();
       }
     mean /= N2up;
     return mean;
   }
}

double Y1up::expected_value()
/* calculates the expected value of the price at current date + delta */
{
  double mean = 0.;
  for (int i = 0; i < NK; i++)
    {
      /* saves */
      int kref = k;
      for(int i = 0; i < k; i++) Lmem1[i] = L[i];
      Scheme();
      mean += e_v_short_traj1();
      /*  loads */
      k = kref;
      for(int i = 0; i < k; i++) L[i] = Lmem1[i];
    }
  mean /= NK;
  return mean;
}

double Y1up::traj1up()
/* one path simulation for Yup^up and Yup_low */
{
  double payo = payoff;
  double e_v1, e_v2, e_v;  
  double actual1 = payo;// for Yup^up
  double actual2 = payo;// for Yup_low
  up = actual1;
  low = actual2;
  for (int i = 0; i < kinit - 1; i++)
    {
      e_v1 = expected_value();// for Yup^up
      e_v2 = expected_value();// for Yup_low
      actual1 += e_v1 - payo;
      actual2 += e_v2 - payo;
      Scheme();
      payo = payoff;
      e_v = e_v_short_traj1();
      actual1 += payo - e_v;
      actual2 += payo - e_v;      
      if (actual1 > up) 
	{
	  up = actual1;
	  low = actual2;
	}
    }
}

double Y1up::alpha_determination()
/* heuristic estimation of alpha */
{
  double meanup, meanlow;
  int N = 2;
  double cl = 0., cu = 0.; 
  for(int j = 1; j <= 5; j++)
    {  
      N *= 2;
      meanup = 0.;
      meanlow = 0.;
      for(int i = 0; i < N; i++)
	{
	  InitialCond(); 
	  traj1up();
	  meanup += up; meanlow += low;  
	} 
      meanup /= N; meanlow /= N;
      cu += meanup; cl += meanlow;
    }
  double alpha = cl/(cu+cl);
  return alpha;
}

void Y1up::estimate()
/* provides the arithmetic mean, standard deviation and a 95% confidence interval for the price upper estimator and the values of alpha for the combined estimator */
{
  double meanup = 0.; 
  double meanlow = 0.; 
  double varup = 0.; 
  double varlow = 0.; 
  for(int i = 0; i < N1up; i++)
    {
      InitialCond(); 
      traj1up();
      meanup += up; 
      meanlow += low;      
      varup += up*up;     
      varlow += low*low;
    }
  meanup /= N1up; varup /= N1up; varup -= meanup*meanup;
  meanlow /= N1up; varlow /= N1up; varlow -= meanlow*meanlow;
  double SDup = sqrt(varup);
  double interup = 1.96*SDup/sqrt(N1up); 
  interup *= 10000; meanup *= 10000;  
  cout << "price upper upper estimator : " << meanup << "  (95 % confidence interval : [" << (meanup - interup) << "," << (meanup + interup) << "]" << " )" << endl;
  double SDlow = sqrt(varlow);
  double interlow = 1.96*SDlow/sqrt(N1up); 
  interlow *= 10000; meanlow *= 10000;
  cout << "price upper lower estimator : " << meanlow << "  (95 % confidence interval : [" << (meanlow - interlow) << "," << (meanlow + interlow) << "]" << " )" << endl;
  double alpha = alpha_determination(); 
  cout << "price upper combined estimator (alpha = " << alpha << "): " << alpha*meanup+(1-alpha)*meanlow << endl;
}


double YAup::short_trajA()
/* simulates Z_{tau_A} */
{ 
      /* saves */
      int kref = k;
      for(int i = 0; i < k; i++) Lmem2[i] = L[i];
      bool b = true; // false : exercice
      for(int i = 0; (i < kref -1 ) && b; i++)
	{     
	  Scheme(); 
	  if (payoff != 0.) // no need to compare  the payoff in this case (frequent in OTM case)
	    {
	      b = (payoff <= H[i]);
	    }
	}     
      /*  loads */
      k = kref;
      for(int i = 0; i < k; i++) L[i] = Lmem2[i];
      return payoff;
}

double YAup::e_v_short_trajA()
/* calculates an estimation of the price with Andersen's algorithm */
{
  if (payoff > H[kinit - k]) 
    {
      return payoff;
    }
  else
    {
      double mean = 0.;
      for (int i = 0; i < N2up; i++)
	{
	  mean += short_trajA();
	}
      mean /= N2up;
      return mean;
    }
}

double YAup::expected_value()
/* calculates the expected value of the price at current date + delta */
{
  double mean = 0.;
  for (int i = 0; i < NK; i++)
    {
      /* saves */
      int kref = k;
      for(int i = 0; i < k; i++) Lmem1[i] = L[i];
      Scheme();
      mean += e_v_short_trajA();
      /*  loads */
      k = kref;
      for(int i = 0; i < k; i++) L[i] = Lmem1[i];
    }
  mean /= NK;
  return mean;
}

void YAup::trajAup()
/* a one path simulation for Yup^up and Yup_low */
{
  double payo = payoff;
  double e_v1, e_v2, e_v;  
  double actual1 = payo; // for Yup^up
  double actual2 = payo; // for Yup_low 
  up = actual1;
  low = actual2;
  for (int i = 0; i < kinit - 1; i++)
    {
      e_v1 = expected_value(); // for Yup^up
      e_v2 = expected_value(); // for Yup_low 
      actual1 += e_v1 - payo;
      actual2 += e_v2 - payo;
      Scheme();
      payo = payoff;
      e_v = e_v_short_trajA();
      actual1 += payo - e_v;
      actual2 += payo - e_v;      
      if (actual1 > up) 
	{
	  up = actual1;
	  low = actual2;
	}
    }
}

double YAup::alpha_determination()
/* heuristic estimation of alpha */
{
  double meanup, meanlow;
  int N = 2;
  double cl = 0., cu = 0.; 
  for(int j = 1; j <= 5; j++)
    {  
      N *= 2;
      meanup = 0.;
      meanlow = 0.;
      for(int i = 0; i < N; i++)
	{
	  InitialCond(); 
	  trajAup();
	  meanup += up; meanlow += low;  
	} 
      cu += meanup; cl += meanlow;
    }
  double alpha = cl/(cu+cl);
  return alpha;
}

void YAup::estimate2()
/* provides the arithmetic mean, standard deviation and a 95% confidence interval for the price upper estimator and the values of alpha for the combined estimator */
{
  Find_H();
  double meanup = 0.;
  double meanlow = 0.; 
  double varup = 0.;
  double varlow = 0.; 
  for(int i = 0; i < N1up; i++)
    {
      InitialCond(); 
      trajAup();
      meanup += up; 
      meanlow += low;      
      varup += up*up;     
      varlow += low*low;
    }
  meanup /= N1up; varup /= N1up; varup -= meanup*meanup;
  meanlow /= N1up; varlow /= N1up; varlow -= meanlow*meanlow;
  double SDup = sqrt(varup);
  double interup = 1.96*SDup/sqrt(N1up);
  interup *= 10000; meanup *= 10000;
  cout << "price upper upper estimator : " << meanup << "  (95 % confidence interval : [" << (meanup - interup) << "," << (meanup + interup) << "]" << " )" << endl;
  double SDlow = sqrt(varlow);
  double interlow = 1.96*SDlow/sqrt(N1up); 
  interlow *= 10000; meanlow *= 10000;
  cout << "price upper lower estimator : " << meanlow << "  (95 % confidence interval : [" << (meanlow - interlow) << "," << (meanlow + interlow) << "]" << " )" << endl;
  double alpha = alpha_determination(); 
  cout << "price upper combined estimator (alpha = " << alpha << "): " << alpha*meanup+(1-alpha)*meanlow << endl;
}
