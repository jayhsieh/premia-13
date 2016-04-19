
/*------------------------------------*/
/*   BERMUDAN SWAPTION PRICER         */
/*    Andersen Algorithm              */
/*                                    */
/*------------------------------------*/
/*  Sonke Blunck, Premia 2005         */
/*------------------------------------*/

#ifndef LMM_BERMUDAPRICE_ANDERSEN_H
#define LMM_BERMUDAPRICE_ANDERSEN_H

#include <valarray>
#include "mathsb.h"


class LMM_1F
// one-factor LIBOR Market Model
{
  int _N;                      // number of LIBORs
  float _delta;                // accrual period = MC time step 
  double _sqrt_delta;          // sqrt of _delta 
  int _eta;                    // current time = (_eta-1)*_delta
  double _PVBP;                // Principal Value of Basis Point
  double _swaprate;            

  std::valarray<double> _L0;          // for the LIBORs at time 0 
  std::valarray<double> _Lt;          // for the LIBORs at time min(t,i*_delta) 
  std::valarray<double> _St;          // for the log of the LIBORs
  std::valarray<float>  _lambda;      // forward volatilities
  std::valarray<double> _sqr_lambda;  
  std::valarray<double> _ZCB;         // for the ZCB price with maturity i*_delta

 public: 
  LMM_1F( int N=20, float delta=0.5, float lambda=0.2, float L0=0.06 );

  std::valarray<double>&  Get_Lt();

  void InitialCond();

  void Scheme( int m ); 
  // simulation from time (_eta-1)*_delta to time _eta*_delta
  // here we define a first-order log-Euler scheme under the spot dynamics 

  void Set_ZCB( int p, int s, int e ); 
  // sets _ZCB[i]=P(p*_delta,i*delta) for i=p,...,e
  // then sets _PVBP and _swaprate (for time p*_delta) according to: 
  // swap start-date = s*_delta   and   swap end-date = e*_delta 

  double SpotNumeraire( int p, std::valarray<double> &L );
  // returns the value of the spot numeraire at time p*_delta 

  double SwaptionCF( int p, std::valarray<double> &L, int s, int e, float K);
  // CF approx. of the value at time t=p*_delta (provided _Lt=L) of the 
  // Europ. swaption with 
  // swap start-date = s*_delta   and   swap end-date = e*_delta 
  // The maturity of the swaption is s*_delta 

  double SwaptionMC( int s, int e, float K, int M );   
  // MC approx. of the present value of the Europ. swaption with 
  // swap start-date = s*_delta   and   swap end-date = e*_delta 
  // The maturity of the swaption is s*_delta 

}; // end of the class LMM_1F 










class FctToMaxim : public NumFct1D
// for the fcts to be maximized in the Andersen algorithm 
{
  int _i;   // index of exercise date at which we want to maximize 
            // the expected discounted value of the Bermudan swaption 
  int _L;   // number of MC paths 

  // all other data members are as in the class Andersen 

  int _s;   // first swaption ex. date = swap start-date = _T[_s]
  int _l;   // last swaption exercise date = _T[_l]
  int _e;   // swap end-date = T[_e]; n° of LIBORs to consider = _e
  float _K; // strike 
  
  int _n;   // number of exercise dates 
  int _m;   // number of swap payment dates 

  std::valarray<double> _H;   // param. of the exercise strategy at T[i] 
  int _strategy;         // Andersen strategy I or II (see paper)

  std::valarray<double> _x; // MC paths: _x[ l*_e*_n + k*_n + j] is the l-th
                       // sample of the LIBOR n° k at _T[j] 

  LMM_1F _LMM;


 public:
  FctToMaxim( int i, int L, int s, int l, int e, float K, 
              LMM_1F LMM, std::valarray<double> H,  int stra,  
	      std::valarray<double> &x ); 

  virtual ~FctToMaxim() {}

  double f( int i, int j, std::valarray<double> &L );
  // CF for the discounted value at time t of the i-th Europ. 
  // component option

  int tilde_b( int i, double H, std::valarray<double> &v );
  // CF for the indicator function of the event "exercising at 
  // _T[i] is optimal" in terms of the discounted values at _T[i] 
  // of the Europ. component options
  // We use the exercise strategy II in the Andersen-paper 

  int b( int i, double H, std::valarray<double> &L );
  // tilde_b written via f as a fct of the driving process x   

  int g( double H, std::valarray<double> &x );
  // CF for the optimal stopping time at _T[_i] 
  
  double F( double H, std::valarray<double> &x );
  // CF (apart from conditioning on _T[_i] !!) for the discounted 
  // value at _T[_i] of the Bermudan option 
  
  double Eval( double H );
  // calculates via MC the expected discounted value at _T[_i] of 
  // the Bermudan swaption 
  // this is the fct which has to be maximized; we multiplied it 
  // with -1 since our Golden Section algo performs a minimization  

};  // end of the class FctToMaxim 







class Andersen 
{
  int _s;   // first swaption ex. date = swap start-date = _T[_s]
  int _l;   // last swaption exercise date = _T[_l]
  int _e;   // swap end-date = T[_e]; n° of LIBORs to consider = _e
  float _K; // strike 
  
  int _n;   // number of exercise dates 
  int _m;   // number of swap payment dates 

  std::valarray<float>  _T;   // swaption exercise dates 
  std::valarray<double> _H;   // param. of the exercise strategy at T[i] 
  int _strategy;         // Andersen strategy I or II (see paper)

  std::valarray<double> _x; // MC paths: _x[ l*_e*_n + k*_n + j ] is the l-th
                       // sample of the LIBOR n° k at _T[j] 

  LMM_1F _LMM;


 public:
  Andersen( int s, int l, int e, float K, float delta, int stra, 
	    float vol, float L0 );

  double f( int i, int j, std::valarray<double> &L );
  // CF for the discounted value at time T[j] of the i-th Europ. 
  // component option

  int tilde_b( int i, std::valarray<double> &v );
  // CF for the indicator function of the event "exercising at 
  // _T[i] is optimal" in terms of the discounted values at _T[i] 
  // of the Europ. component options
  // We use the exercise strategy II in the Andersen-paper 

  int b( int i, std::valarray<double> &x );
  // tilde_b written via f as a fct of the driving process x   

  int g( int i, std::valarray<double> &x );
  // CF for the optimal stopping time at _T[i] 

  double F( int i, std::valarray<double> &x );
  // CF (apart from conditioning on _T[i] !!) for the discounted 
  // value at _T[i] of the Bermudan option 
  
  double BermSwaptionMC( int L=10000, int M=50000 );
  // calculates via Andersen-MC the (discounted) value at time 0 
  // of the Bermudan swaption 
  // (value at 0 = discounted value at 0  since SpotNum. at 0 = 1)

};  // end of the class Andersen 


double lmm_swaption_payer_bermudan_andersen_pricer(float tenor, float swaptionMat, float swapMat, float K, float flatInitialValue, float vol, int numberMCPaths1, int numberMCPaths2, int AndersenStrategy );


#endif     
