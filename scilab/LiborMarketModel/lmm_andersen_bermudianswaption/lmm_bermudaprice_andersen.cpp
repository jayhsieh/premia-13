/*------------------------------------*/
/*   BERMUDAN SWAPTION PRICER         */
/*    Andersen Algorithm              */
/*                                    */
/*------------------------------------*/
/*  Sonke Blunck, Premia 2005         */
/*------------------------------------*/

#include <iostream>
#include <cmath>
#include "mathsb.h"


#include "lmm_bermudaprice_andersen.h"

using namespace std;



static double TWO_PI = 6.28318530717958623; 
static double SQRT_TWO_PI_INV = 1/sqrt(TWO_PI);  



double U()
{
  return drand48(); 
}



double G()
{
  static double sqrt_term, u;
  static int Sin=1; 

  if (Sin) 
    {
      sqrt_term = sqrt( (-2.)*log(U()) ); 
      u = U();
      Sin=0;
      return sqrt_term * sin(TWO_PI*u);
    }
  else
    {
      Sin=1;
      return sqrt_term * cos(TWO_PI*u);
    }
}



int intapprox(double t)
{
  int n=(int)t;

  if (t-n<=0.5) return n; else return n+1;
}






//////////////////////////////////////////////////
//
// class LMM_1F  (one-factor LIBOR Market Model) 
//
//////////////////////////////////////////////////

LMM_1F::LMM_1F( int N, float delta, float vol, float L0 ):
    _N(N), _delta(delta)
  {
    _sqrt_delta = sqrt(_delta);

    _lambda.resize(_N); _lambda=vol;  
    _sqr_lambda.resize(_N); _sqr_lambda=SQR(vol);
    _L0.resize(_N); _L0=L0; 
    _Lt.resize(_N);
    _St.resize(_N);
    _ZCB.resize(_N+1);

  } // end of the constructor



valarray<double>& LMM_1F::Get_Lt()
{ return _Lt; }



void LMM_1F::InitialCond()
  { 
    _Lt=_L0; 
    _St=log(_L0); 
    _eta=1; 
  }




void LMM_1F::Scheme( int m )
  // simulation from time (_eta-1)*_delta to time _eta*_delta
  // here we define a first-order log-Euler scheme under the spot dynamics 
  {
    double drift;
    double deltaW = _sqrt_delta*G();

    for (int k=_eta; k<m; k++)
      {
	drift=0.; 
	for (int j=_eta; j<=k; j++)  
	  drift+=_Lt[j]*_lambda[j]/(1+_delta*_Lt[j]); 
	drift=_delta*_lambda[k]*drift - _sqr_lambda[k]/2.;

	_St[k]+=drift*_delta + _lambda[k]*deltaW;
      }

    _Lt=exp(_St);
    _eta++;
  }



void LMM_1F::Set_ZCB( int p, int s, int e )
  // sets _ZCB[i]=P(p*_delta,i*delta) for i=p,...,e
  // then sets _PVBP and _swaprate (for time p*_delta) according to: 
  // swap start-date = s*_delta   and   swap end-date = e*_delta 
  {
    _ZCB[p]=1.; 
    for (int i=p+1; i<=e; i++)
      _ZCB[i] = _ZCB[i-1] / (1+_delta*_Lt[i-1]); 

    _PVBP=0.;
    for (int i=s+1; i<=e; i++)  _PVBP+=_ZCB[i];
    _PVBP*=_delta;

    _swaprate=(_ZCB[s]-_ZCB[e])/_PVBP;
  }



double LMM_1F::SpotNumeraire( int p, valarray<double> &L )
  // returns the value of the spot numeraire at time p*_delta 
  {
    double prod=1.;

    for (int j=0; j<p; j++)  prod*=1+_delta*L[j];
    return prod;
  }



double LMM_1F::SwaptionCF( int p, valarray<double> &L, int s,
			   int e, float K)
  // CF approx. of the value at time t=p*_delta (provided _Lt=L) of the 
  // Europ. swaption with 
  // swap start-date = s*_delta   and   swap end-date = e*_delta 
  // The maturity of the swaption is s*_delta 
  {
    double vol=0., Psum=0., Pterm;

    _Lt=L;
    Set_ZCB( p, s, e ); 

    // if p=s (i.e. we are at maturity !) return the payoff 
    if (p==s)  return _PVBP * MAX( _swaprate-K , 0. ); 


    Pterm = _ZCB[e] / (_ZCB[s] - _ZCB[e]); 

    for (int j=e-1; j>=s; j--)
      {
	Psum+=_ZCB[j+1];
	vol+=_lambda[j]*L[j]/(1.+_delta*L[j])*(Pterm + _delta*Psum/_PVBP);
      }
    vol*=_delta;


    double sigma, logterm, dplus, dminus; // param. for BS-formula
    sigma = sqrt(s*_delta - p*_delta)*vol;
    logterm = log(_swaprate/K); 
    dplus  = (logterm + SQR(sigma)/2.) / sigma; 
    dminus = (logterm - SQR(sigma)/2.) / sigma;

    return _PVBP * ( _swaprate*Cumul_Normal(dplus) - 
		           K*Cumul_Normal(dminus) );
  }


 
double LMM_1F::SwaptionMC( int s, int e, float K, int M )   
  // MC approx. of the present value of the Europ. swaption with 
  // swap start-date = s*_delta   and   swap end-date = e*_delta 
  // The maturity of the swaption is s*_delta 
  {
    int l,j;
    double disc_payoff, var_estim, sum=0., sqr_sum=0.;

    for (l=0; l<M; l++)
      {
	InitialCond();
	for (j=0; j<s; j++)  Scheme(e);

	Set_ZCB( s, s, e );
	disc_payoff = _PVBP / SpotNumeraire(s,_Lt) * MAX( _swaprate-K , 0. );
	sum+=disc_payoff;
	sqr_sum+=SQR(disc_payoff);
      }

    var_estim = (sqr_sum - SQR(sum)/(double)M) / (double)(M-1); 
    // cout << endl << "95% conf. interval = " 
    //      << 10000*1.96*sqrt(var_estim/M) << endl; 

    return sum/(double)M;
  }





//////////////////////////////////////////////////
//
// class FctToMaxim 
//
//////////////////////////////////////////////////

FctToMaxim::FctToMaxim( int i, int L, int s, int l, int e, float K, 
			LMM_1F LMM, valarray<double> H,  int stra,
			valarray<double> &x ):
    NumFct1D(),
    _i(i), _L(L), _s(s), _l(l), _e(e), _K(K), _H(H), _strategy(stra)
    {  
      _LMM=LMM;
      _n=_l-_s+1;
      _m=_e-_s;
      
      _x.resize(L*_n*_e); 
      for (int j=0;  j<L*_n*_e; j++)  _x[j]=x[j]; 

    }  // end of the constructor



double FctToMaxim::f( int i, int j, valarray<double> &L )
  // CF for the discounted value at time t of the i-th Europ. 
  // component option
  { return _LMM.SwaptionCF(j+_s,L,i+_s,_e,_K) / _LMM.SpotNumeraire(j+_s,L); }



int FctToMaxim::tilde_b( int i, double H, valarray<double> &v )
  // CF for the indicator function of the event "exercising at 
  // _T[i] is optimal" in terms of the discounted values at _T[i] 
  // of the Europ. component options
  // We use the exercise strategy II in the Andersen-paper 
  {
    if (i<_i)  cout << "Pbm in tilde_b !!" << endl; 

    // case i=_n-1 
    if (i==_n-1)
      {
	if (v[_n-1] <= 0.)  return 0; 
	return 1;
      }
    
    // case i=_i  (here H is used !)
    if (i==_i)
      {
	if (v[i] <= H)  return 0;
	if (_strategy==2)
	  {
	    for (int j=i+1; j<_n; j++)  if (v[i] <= v[j]) return 0;
	  }
	return 1;
      }

    // other cases 
    if (v[i] <= _H[i])  return 0;
    if (_strategy==2)
      {
	for (int j=i+1; j<_n; j++)  if (v[i] <= v[j]) return 0;
      }
    return 1;
    
  }


int FctToMaxim::b( int i, double H, valarray<double> &L )
  // tilde_b written via f as a fct of the driving process x 
  {
    valarray<double> v(0.,_n);

    for (int j=i; j<_n; j++)  v[j] = f(j,i,L);
    return tilde_b(i,H,v);
  } 
  

int FctToMaxim::g( double H, valarray<double> &x )
// CF for the optimal stopping time at _T[_i] 
// x[k*_n+j] = LIBOR n° k at time _T[j] 
{ 
  valarray<double> L(0.,_e);

  for (int j=_i; j<_n; j++)
    {
      for (int k=0; k<_e; k++)  L[k]=x[k*_n+j];
      if (b(j,H,L)) return j;
    }
  return _n;  // no-exercise case !!    
}

  
double FctToMaxim::F( double H, valarray<double> &x )
// CF (apart from conditioning on _T[_i] !!) for the discounted 
// value at _T[_i] of the Bermudan option 
  {
    valarray<double> L(0.,_e); 
    int tau = g(H,x); 

    if (tau >= _n)  return 0.;  // no-exercise case !! 
        
    for (int k=0; k<_e; k++)  L[k]=x[k*_n+tau];
    return f(tau, tau, L);
  }
  

  
double FctToMaxim::Eval( double H )
// calculates via MC the expected discounted value at _T[_i] of 
// the Bermudan swaption 
// this is the fct which has to be maximized; we multiplied it 
// with -1 since our Golden Section algo performs a minimization  
  {
    double sum=0.; 
    int j,k,l;
    valarray<double> x(0.,_e*_n); 
  
    for (l=0; l<_L; l++)
      {
	for (k=0; k<_e; k++)
	  for (j=0; j<_n; j++)
	    x[k*_n+j] = _x[l*_e*_n + k*_n + j]; // l-th sample at T[j]
	sum+=F(H,x);
      }

    //    cout << "H=" << H << "   Eval(H)=" << sum << endl; 

    return -sum; 
  }







//////////////////////////////////////////////////
//
// class Andersen 
//
//////////////////////////////////////////////////

Andersen::Andersen( int s, int l, int e, float K, float delta, 
		    int stra, float vol, float L0 ):
  _s(s), _l(l), _e(e), _K(K), _strategy(stra), _LMM(e,delta,vol,L0)
  {
    _n=_l-_s+1; 
    _m=_e-_s;
    _H.resize(_n-1);
    _T.resize(_n); for (int i=0; i<_n; i++)  _T[i]=(i+_s)*delta; 
    _x.resize(1);

  } // end of the constructor



double Andersen::f( int i, int j, valarray<double> &L )
  // CF for the discounted value at time T[j] of the i-th Europ. 
  // component option
  { return _LMM.SwaptionCF(j+_s,L,i+_s,_e,_K) / _LMM.SpotNumeraire(j+_s,L); }
  


int Andersen::tilde_b( int i, valarray<double> &v )
  // CF for the indicator function of the event "exercising at 
  // _T[i] is optimal" in terms of the discounted values at _T[i] 
  // of the Europ. component options
  // We use the exercise strategy II in the Andersen-paper 
  {
    // case i=_n-1 
    if (i==_n-1)
      {
	if (v[_n-1]>0.)  return 1; else return 0;
      }

    // other cases 
    if (v[i] <= _H[i])  return 0;
    if (_strategy==2)
      { 
	for (int j=i+1; j<_n; j++)  if (v[i] <= v[j]) return 0; 	
      }
    return 1;
  }



int Andersen::b( int i, valarray<double> &L )
// tilde_b written via f as a fct of the driving process x 
  {
    valarray<double> v(0.,_n);

    for (int j=i; j<_n; j++)  v[j] = f(j,i,L);
    return tilde_b(i,v);
  } 

  

int Andersen::g( int i, valarray<double> &x )
// CF for the optimal stopping time at _T[_i] 
// x[k*_n+j] = LIBOR n° k at time _T[j] 
{ 
  valarray<double> L(0.,_e);

  for (int j=i; j<_n; j++)
    {
      for (int k=0; k<_e; k++)  L[k]=x[k*_n+j];
      if (b(j,L)) return j;
    }
  return _n;  // no-exercise case !!    
}



double Andersen::F( int i, valarray<double> &x )
// CF (apart from conditioning on _T[i] !!) for the discounted 
// value at _T[i] of the Bermudan option 
{
  valarray<double> L(0.,_e);
  int tau = g(i,x);

  if (tau >= _n)  return 0.;  // no-exercise case !! 
        
  for (int k=0; k<_e; k++)  L[k]=x[k*_n+tau];
  return f(tau, tau, L);
}



double Andersen::BermSwaptionMC( int L, int M )
// calculates via Andersen-MC the (discounted) value at time 0 
// of the Bermudan swaption 
// (value at 0 = discounted value at 0  since SpotNum. at 0 = 1)
{
  int i,j,k,l,m;
  valarray<double> Libor(0.,_e);

  // Construction of L MC paths 
  _x.resize(L*_n*_e); 

  for (l=0; l<L; l++)
    {
      _LMM.InitialCond();
      for (j=0; j<_s-1; j++)  _LMM.Scheme(_e);;

      for (j=0; j<_n; j++)
	{
	  _LMM.Scheme(_e);  // now we are at time T[j] 
	  Libor = _LMM.Get_Lt();
	  for (k=0; k<_e; k++)
	    _x[l*_e*_n + k*_n + j] = Libor[k];  
          // l-th sample at time T[j]
	}
    }
 
  // Computation of the _H[i] via optimization 
  FctToMaxim A(_n-2,L,_s,_l,_e,_K,_LMM,_H,_strategy,_x);
  double ax,bx; // start values for Golden Section Search 
  _H=0.; 
  for (i=_n-2; i>=0; i--)
    {
      A=FctToMaxim(i,L,_s,_l,_e,_K,_LMM,_H,_strategy,_x); 
      ax=_H[ MIN(i+1,_n-2) ]; bx=ax+0.001; 
      GoldenSectionMin1D( A, ax, bx, _H[i] );
      // cout << "_H[" << i << "]=" << _H[i] << endl;
    }
  _x.resize(0);


  // Computation of the current discounted price of the Berm. Swaption 
  // via M new MC paths 
  double disc_payoff, var_estim, sum=0., sqr_sum=0.;
  valarray<double> x(0.,_e*_n);

  for (m=0; m<M; m++)
    {
      _LMM.InitialCond();
      for (j=0; j<_s-1; j++)  _LMM.Scheme(_e);

      for (j=0; j<_n; j++)
	{
	  _LMM.Scheme(_e);  // now we are at time T[j] 
	  Libor = _LMM.Get_Lt();
	  for (k=0; k<_e; k++)  x[k*_n+j]=Libor[k];
	}

      disc_payoff=F(0,x); 
      sum+=disc_payoff;
      sqr_sum+=SQR(disc_payoff);
    }

  //cout << "price = " << 10000*sum/(double)M << endl;  

  var_estim = (sqr_sum - SQR(sum)/(double)M) / (double)(M-1);
  // cout << endl << "95% conf. interval = " 
  //      << 10000*1.96*sqrt(var_estim/M) << endl; 

  return sum/(double)M;

}





double lmm_swaption_payer_bermudan_andersen_pricer(float tenor, float swaptionMat, float swapMat, float K, float flatInitialValue, float vol, int numberMCPaths1, int numberMCPaths2, int AndersenStrategy )
{
  int s=intapprox( swaptionMat/tenor );  // Start date of swaption and swap
  int e=intapprox( swapMat/tenor );      // End date of swap
  int l=e-1;                             // Last swaption ex. date 

  Andersen Ander(s,l,e,K,tenor,AndersenStrategy,vol,flatInitialValue);

  return Ander.BermSwaptionMC(numberMCPaths1, numberMCPaths2);
}






// int main()
// {

//   double p;
//   float tenor=0.5;              //usually 3 months or 6 months 
//   float swaptionMat=1.0;        //(years)
//   float swapMat=4.0;            //(years)
//   float K=0.06;                 //strike
//   float flatInitialValue=0.06;  //for the SDE of the LIBORs
//   float vol=0.2;                //for the SDE of the LIBORs
//   int numberMCPaths1=10000;     //for the comput. of the parameters
//   int numberMCPaths2=50000;     //for the comput. of the price
//   int AndersenStrategy=1;       //must be 1 or 2 See Documentation and Paper

//   srand48(1);

//   p=lmm_swaption_payer_bermudan_andersen_pricer(tenor, swaptionMat, swapMat, K, flatInitialValue, vol, numberMCPaths1, numberMCPaths2, AndersenStrategy); 

//   printf("Bermudan swaption with Andersen algorithm is %f\n",p*10000);
 
//   return(EXIT_SUCCESS);
// }



