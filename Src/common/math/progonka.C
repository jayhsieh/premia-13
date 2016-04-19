#include "progonka.h"
#include "numerics.h"

vector<double> progonka(int N, const double low, 
			const double diag, const double up, 
			const vector<double> & right_side)
{
  vector<double> A(N-1), B(N-1), result(N);
  
  // descente 
  A[0] = -up/diag;
  B[0] = right_side[0]/diag;
  double denom;
  for(int i=1;i<N-1;i++)
    { 
      denom = diag + low*A[i-1];  
      A[i] = -up/denom;
      B[i] = (right_side[i]-low*B[i-1])/denom;
    }
  
  // remontee
  result[N-1] = (right_side[N-1] - low*B[N-2])/(diag + low*A[N-2]);
  for(int i=N-2;i>=0;i--)
    {
      result[i] = A[i]*result[i+1] + B[i];
    }
  
  return result;
}

        
        
        
        
        
        
        
