#include "stack-c.h"

extern double  implied_volatility(double price,double r, double S0, double T, double K,double error);


int implied(char *fname)
{
  int m1,n1,l1,m2,n2,l2,m3,n3,l3,m4,n4,l4,m5,n5,l5,n6,l6,m6;
  int minlhs=1,maxlhs=1,minrhs=6,maxrhs=6;
  double price,r,S0,T,K,error;
  int m,n,y;
  

  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);

  GetRhsVar(1,"d",&m1,&n1,&l1);
  GetRhsVar(2,"d",&m2,&n2,&l2);
  GetRhsVar(3,"d",&m3,&n3,&l3);
  GetRhsVar(4,"d",&m4,&n4,&l4);
  GetRhsVar(5,"d",&m5,&n5,&l5);
  GetRhsVar(6,"d",&m6,&n6,&l6);

  if (m1 * n1 !=1 || m2 * n2 !=1 || m3 * n3 !=1 || m4 * n4 !=1 || m5 * n5 !=1 || m6 * n6 !=1){
	cerro("Les arguments de Implied doivent etre scalaires");
	return 0;
  }
  
  price=*stk(l1);
  r=*stk(l2);
  S0=*stk(l3);
  T=*stk(l4);
  K=*stk(l5);
  error=*stk(l6);
  m=1;
  n=1;
  
  CreateVar(7,"d",&m,&n,&y);
  *stk(y) = implied_volatility(price,r,S0,T,K,error);
  
  LhsVar(1)=7;
  return 0;
}

