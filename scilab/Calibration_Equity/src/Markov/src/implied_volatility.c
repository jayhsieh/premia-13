#include <math.h>

double LoiNormale(double x)
{
  const double p= 0.2316419;
  const double b1= 0.319381530;
  const double b2= -0.356563782;
  const double b3= 1.781477937;
  const double b4= -1.821255978;
  const double b5= 1.330274429;
  const double one_over_twopi= 0.39894228;

  double t;

  if(x >= 0.0)
  {
    t = 1.0 / ( 1.0 + p * x );
    return (1.0 - one_over_twopi*exp(-x*x/2.)*t*(t*(t*(t*(t*b5+b4)+b3)+b2)+b1));
  } 
  else 
  {/* x < 0 */
    t = 1.0 / ( 1.0 - p * x );
    return (one_over_twopi*exp(-x*x/2.)*t*(t*(t*(t*(t*b5+b4)+b3)+b2)+b1));
  }
}

double implied_volatility(double price,double r, double S0, double T, double K,double error)
{

  int  lc_step, ln_loop, pi;
  double vol, dv;
  double d1, d2;
  double price_error, risk_lessT, vega;
  double real_price;
  double dv_old, dv_best, vol_best;  // //

  double PI=3.141592;
  double ome=1.0;
  double mu=0.0;
  
  lc_step = 0;
  // init parameters
  // start point of the algo with a given volatility
  vol =  1.0 ; // //   0.300;
  dv  = error + 1;
  pi = PI;
  ln_loop = 0;
  dv_best = 1000;
  vol_best = vol;
  real_price = price; // market price of the option

  lc_step = 10;
  // loop
  while ( fabs(dv) > error)
	{
	  lc_step = 15;
	  // stopping point
	  dv_old = dv;
	  // phases d1 and and d2
	  d1 = ome * ( log(S0 / K ) +
				   (mu + 0.5*vol*vol )  * T
				   );
	  d1 = d1 / (vol * sqrt(T));
	  d2 = d1 - ome * vol * sqrt(T);

	  // price calculation with volatility vol
	  lc_step = 20;
	  price_error =  S0 * exp (mu * T)  *  LoiNormale(d1)
		- K * LoiNormale(d2);
	  risk_lessT = exp(   - T *r);
	  price_error = ome *  risk_lessT * price_error -
		price;
	  // vega calculation
	  lc_step = 30;
	  vega = risk_lessT * S0 * exp (mu * T) * sqrt(T)
		* (1/(sqrt(2*pi))) * exp (- 0.50 * d1*d1);

	  // correction
	  lc_step = 40;
	  dv  = ( price_error / vega);

	  //  if vega becomes too low then slowdown the decrease rate
	  if ((  fabs(dv / dv_old ) < 1 ))
		{
		  lc_step = 50;
		  // better than best dv ?
		  if  ((fabs(dv) > fabs(dv_best))  && (fabs(vega) < 1))
			{
			  lc_step = 55;
			  //  stop research
			  vol = vol_best;
			  dv = 0;
			}
		  else
			{
			  lc_step = 60;
			  // research with no correction
			  vol = vol - dv;
			  if ( fabs(dv) < fabs(dv_best) )
				{
				  // improving record
				  dv_best = dv;
				  vol_best = vol;
				}
			}
		  lc_step = 70;
		}
	  else
		{
		  // correction
		  lc_step = 80;
		  vol = vol + dv_old - dv_old*dv_old;
		  dv = dv_old;
		  if (vol >= 2.5)
			{ vol = vol_best;
			dv = vol_best - dv_best*dv_best;}
		  //  vol_best = 0;
		  //  dv = 0;
		}

	  lc_step = 90;
	  ln_loop = ln_loop + 1;
	  // too much loop > stop
	  if (ln_loop > 1000)
		{ dv = 0;}
	}

  //stat
  // //
  return vol_best;
}


