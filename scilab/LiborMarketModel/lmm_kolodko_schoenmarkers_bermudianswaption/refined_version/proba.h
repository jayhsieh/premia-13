#define TWO_PI 6.2831853

double U()
{
  return drand48(); 
}

/* Normal distribution */
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

/*One-Dimensional Normal Law. Cumulative distribution function. */
/*Abramowitz, Milton et Stegun, Handbook of MathematicalFunctions, 1968, Dover Publication, New York, page 932 (26.2.18).Precision 10-7*/
double Cumul_Normal(double x)
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
      return (1.0 - one_over_twopi * exp( -x * x / 2.0 ) * t * ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    } 
  else 
    {/* x < 0 */
      t = 1.0 / ( 1.0 - p * x );
      return ( one_over_twopi * exp( -x * x / 2.0 ) * t * ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}

