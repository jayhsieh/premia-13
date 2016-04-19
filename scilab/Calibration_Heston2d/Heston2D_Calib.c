
/****************************************************************
The goal:
---------
Knowing all the marginal (1D) parameters of a 2D Heston model,
this routine fixes the value of the stocks' correlation from
the price of a spread option.

The code:
---------
This C code is parallelized using OpenMPs pragmas in order to 
make it more efficient.

The author:
-----------
This code was written by: Lokman A. Abbas Turki
It is included to the Premia library

If any bug was detected, you can contact the author by sending
an email to:

lokman.abbasturki@gmail.com

*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define SQR(x) ((x) * (x))
#define MAX(A,B) ( (A) > (B) ? (A):(B) )
#define MIN(A,B) ( (A) < (B) ? (A):(B) )
#define ABS(x) ( (x) > (0) ? (x):(-x) )
#define PI 3.14159265358979

///////////////////////////////////////////////////////////////////////////////
// L'Eucuyer CMRG Matrix Values for the Random Numbers Generator
///////////////////////////////////////////////////////////////////////////////
// First MRG 
#define a12 63308
#define a13 -183326
#define q12 33921
#define q13 11714
#define r12 12979
#define r13 2883

// Second MRG 
#define a21 86098
#define a23 -539608
#define q21 24919
#define q23 3976
#define r21 7417
#define r23 2071

// Normalization variables
#define Invmp 4.656612873077393e-10
#define two17   131072.0
#define two53   9007199254740992.0

// Np parameters
#define  pNp  0.2316419
#define  b1 0.319381530
#define  b2 -0.356563782
#define  b3 1.781477937
#define  b4 -1.821255978
#define  b5 1.330274429
#define  one_over_twopi 0.39894228   

// Matrices needed for initialization of the RNG
double A1[3][3], A2[3][3];

// Cholesky decomposition of the correlation matrix
float L[4][4];

///////////////////////////////////////////////////////////////////////////////
// Dynamic variables
///////////////////////////////////////////////////////////////////////////////
// The computation parameters 
float *TPE;
float **Spot1, **Var1;
float **Spot2, **Var2;

// The random numbers 
int **CMRG; 

///////////////////////////////////////////////////////////////////////////////
// Memory allocation of the computation parameters 
///////////////////////////////////////////////////////////////////////////////
static void memory_allocation(int Ntraj, int Nst)
{
  int i;


// The computation parameters 
///////////////////////////////////////////////////////////////////////////////
  // Payoff
  TPE=(float *)malloc((Ntraj)*sizeof(float));

  // First asset
  Spot1=(float **)calloc(Nst,sizeof(float *));
  for (i=0;i<Nst;i++)
    Spot1[i]=(float *)calloc(Ntraj,sizeof(float));
  Var1=(float **)calloc(Nst,sizeof(float *));
  for (i=0;i<Nst;i++)
    Var1[i]=(float *)calloc(Ntraj,sizeof(float));

  // Second asset
  Spot2=(float **)calloc(Nst,sizeof(float *));
  for (i=0;i<Nst;i++)
    Spot2[i]=(float *)calloc(Ntraj,sizeof(float));
  Var2=(float **)calloc(Nst,sizeof(float *));
  for (i=0;i<Nst;i++)
    Var2[i]=(float *)calloc(Ntraj,sizeof(float));


// The random numbers generator parameters
///////////////////////////////////////////////////////////////////////////////
  CMRG=(int **)calloc(6,sizeof(int *));
  for (i=0;i<6;i++)
    CMRG[i]=(int *)calloc(Ntraj,sizeof(int));

   A1[0][0] = 2119409629.0; A1[0][1] =  302707381.0; A1[0][2] =  655487731.0;
   A1[1][0] =  946145520.0; A1[1][1] = 1762689149.0; A1[1][2] =  302707381.0;
   A1[2][0] = 1139076568.0; A1[2][1] =  600956040.0; A1[2][2] = 1762689149.0;

   A2[0][0] = 1705536521.0; A2[0][1] = 1409357255.0; A2[0][2] = 1489714515.0;
   A2[1][0] = 1443451163.0; A2[1][1] = 1705536521.0; A2[1][2] = 1556328147.0;
   A2[2][0] = 1624922073.0; A2[2][1] = 1443451163.0; A2[2][2] =  130172503.0;
}

///////////////////////////////////////////////////////////////////////////////
// Freeing the memory 
///////////////////////////////////////////////////////////////////////////////
static void free_memory(int Nst)
{
  int i;

// The computation parameters 
///////////////////////////////////////////////////////////////////////////////
  // Payoff
  free(TPE);
  // First asset
  for (i=0;i<Nst;i++)
    free(Spot1[i]);
  free(Spot1);
  for (i=0;i<Nst;i++)
    free(Var1[i]);
  free(Var1);

  // Second asset
  for (i=0;i<Nst;i++)
    free(Spot2[i]);
  free(Spot2);
  for (i=0;i<Nst;i++)
    free(Var2[i]);
  free(Var2);

// The random numbers generator parameters
///////////////////////////////////////////////////////////////////////////////
  for (i=0;i<6;i++)
    free(CMRG[i]);
  free(CMRG);

  
}


////////////////////////////////////////////////////////////////////////////
// Cholesky decomposition 4X4 for specific correlatin structure 
////////////////////////////////////////////////////////////////////////////
static void CholDec(float r, float r1, float r2){

	L[0][0]=1.0; 

	L[1][0]=r;
	L[1][1]=sqrt(1.0-r*r);

	L[2][0]=r1;
	L[2][1]=0.0;
	L[2][2]=sqrt(1.0-r1*r1);

	L[3][0]=r*r2;
	L[3][1]=sqrt(1.0-r*r)*r2;
	L[3][2]=0.0;
	L[3][3]=sqrt(1.0-r2*r2);
}



////////////////////////////////////////////////////////////////////////////
// This function generates uniform or gaussian random variable 
// if dist 
////////////////////////////////////////////////////////////////////////////
static double ComputeUnifOrGauss(int dist, int i)
{
 int k;                     // loop index
 int h, p12, p13, p21, p23;    // local variables
 unsigned int inter[2];
 double aleat1, aleat2, loc;

 if(dist!=1){
   if(dist!=2){
	 printf("The value of dist must be 1 or 2");
	 getchar();
   }
 }

 for(k=0;k<dist;k++){
	    	
	const int m1 = 2147483647;
    const int m2 = 2145483479;
		
	// First Component 
	h = CMRG[0][i]/q13; 
	p13 = a13*(h*q13-CMRG[0][i])-h*r13;
		
	h = CMRG[1][i]/q12; 
	p12 = a12*(CMRG[1][i]-h*q12)-h*r12;
		
	if(p13 < 0){
	  p13 = p13 + m1;
	}
	if(p12 < 0){
	  p12 = p12 + m1;
	}
		  
	CMRG[0][i] = CMRG[1][i];
	CMRG[1][i] = CMRG[2][i];
		 
	if((p12 - p13) < 0){
	  CMRG[2][i] = p12 - p13 + m1;  
	}else{CMRG[2][i] = p12 - p13;}
		  
	// Second Component 
	h = CMRG[3][i]/q23; 
	p23 = a23*(h*q23-CMRG[3][i])-h*r23;
		
	h = CMRG[5][i]/q21; 
	p21 = a21*(CMRG[5][i]-h*q21)-h*r21;
		
	if(p23 < 0){
	  p23 = p23 + m2;
	}
	if(p12 < 0){
	  p21 = p21 + m2;
	}
		  
	CMRG[3][i] = CMRG[4][i];
	CMRG[4][i] = CMRG[5][i];
		
	if((p21 - p23) < 0){
	  CMRG[5][i] = p21 - p23 + m2;  
	}else{CMRG[5][i] = p21 - p23;} 
	       
	// Combination
	if(CMRG[2][i] < CMRG[5][i]){
	  inter[k] = CMRG[2][i] - CMRG[5][i] + m1;
	}else {inter[k] = CMRG[2][i] - CMRG[5][i];}    
		   
	if(inter[k]==0){inter[k] = m1;}    
 }

  if(dist==2){
    aleat1 = Invmp*inter[0];
    aleat2 = Invmp*inter[1];
  
    if (aleat1 < 1.45e-6){
	  loc = sqrt(-2*log((double)0.00001))*cos((double)2*PI*aleat2);
	} else {
	     if (aleat1 > 0.99999){
	       loc = 0.0;
		 } else{
	      loc = sqrt(-2*log((double)aleat1))*cos((double)2*PI*aleat2);
		 }
	}
	return loc;
  } else {return Invmp*(double)inter[0];}
  
}




static void HestonSimulation_Modif (float rnr, int Ntraj, float T, int Nst,  
								    float S10, float d1, float v10, float k1, 
								    float o1, float n1, float r1,  
								    float S20, float d2, float v20, float k2, 
								    float o2, float n2, float r2,
								    float r)
{
    long i, m;
    
	float g1, g3, mu1;
	float g2, g4, mu2;

    float X1v, X1s, X2v, X2s;

	float cb1, cb2, cb3, cb4;

	float dt = sqrt((float)T/Nst);

	mu1 = rnr-d1;
	mu2 = rnr-d2;

	CholDec(r, r1, r2);

#pragma omp parallel for private(m,i,X1v,X1s,X2v,X2s,g1,g3,g2,g4,cb1,cb2,cb3,cb4)
    
	for (m=0; m<Ntraj; m++)
    {		
		g1 = ComputeUnifOrGauss(2,m);
		g2 = ComputeUnifOrGauss(2,m);
		g3 = ComputeUnifOrGauss(2,m);
		g4 = ComputeUnifOrGauss(2,m);

		cb1 = g1; 
        cb2 = L[1][0]*g1 + L[1][1]*g2;
		cb3 = L[2][0]*g1 + L[2][1]*g2 + L[2][2]*g3;
		cb4 = L[3][0]*g1 + L[3][1]*g2 + L[3][2]*g3 + L[3][3]*g4;
        
		// First asset
		X1v = SQR(sqrt(v10) + 0.5*n1*dt*cb3) - k1*(v10-o1)*dt*dt -
			  0.25*n1*n1*dt*dt; // Variance
        X1s = (mu1-0.5*v10)*dt*dt + sqrt(v10)*dt*cb1; // Spot

		//Var1[0][m] = X1v;
        if(X1v < 0.0){Var1[0][m] = 0.0;
		}else{Var1[0][m] = X1v;}
		Spot1[0][m] = S10*exp(X1s);

		// Second asset
		X2v = SQR(sqrt(v20) + 0.5*n2*dt*cb4) - k2*(v20-o2)*dt*dt -
			  0.25*n2*n2*dt*dt; // Variance
        X2s = (mu2-0.5*v20)*dt*dt + sqrt(v20)*dt*cb2; // Spot

		//Var2[0][m] = X2v;
        if(X2v < 0.0){Var2[0][m] = 0.0;
		}else{Var2[0][m] = X2v;}
		Spot2[0][m] = S20*exp(X2s);

        for (i=1 ; i<=Nst-1 ; i++)
        {
			g1 = ComputeUnifOrGauss(2,m);
			g2 = ComputeUnifOrGauss(2,m);
			g3 = ComputeUnifOrGauss(2,m);
			g4 = ComputeUnifOrGauss(2,m);

			cb1 = g1; 
			cb2 = L[1][0]*g1 + L[1][1]*g2;
			cb3 = L[2][0]*g1 + L[2][1]*g2 + L[2][2]*g3;
			cb4 = L[3][0]*g1 + L[3][1]*g2 + L[3][2]*g3 + L[3][3]*g4;
	        
			// First asset
			X1v = SQR(sqrt(Var1[i-1][m]) + 0.5*n1*dt*cb3) - 
				  k1*(Var1[i-1][m]-o1)*dt*dt -
				  0.25*n1*n1*dt*dt; // Variance
			X1s = (mu1-0.5*Var1[i-1][m])*dt*dt + 
				  sqrt(Var1[i-1][m])*dt*cb1; // Spot

			//Var1[i][m] = X1v;
			if(X1v < 0.0){Var1[i][m] = 0.0;
			}else{Var1[i][m] = X1v;}
			Spot1[i][m] = Spot1[i-1][m]*exp(X1s);

			// Second asset
			X2v = SQR(sqrt(Var2[i-1][m]) + 0.5*n2*dt*cb4) - 
				  k2*(Var2[i-1][m]-o2)*dt*dt -
				  0.25*n2*n2*dt*dt; // Variance
			X2s = (mu2-0.5*Var2[i-1][m])*dt*dt + 
				  sqrt(Var2[i-1][m])*dt*cb2; // Spot

			//Var2[i][m] = X2v;
			if(X2v < 0.0){Var2[i][m] = 0.0;
			}else{Var2[i][m] = X2v;}
			Spot2[i][m] = Spot2[i-1][m]*exp(X2s);
            
        }

    }

}


 
static void ComputePay(float Str, int Ntraj, int Jindex){

	int i;
	#pragma omp parallel for private(i)
	for (i = 0; i < Ntraj; i++) { 
	   TPE[i] = MAX(Spot1[Jindex][i]-Spot2[Jindex][i]+Str,0.0); 			
	}  
}



static void European_Heston(int Ntraj, int Nst, float rnr, float Str, float T, 
							float d1, float S10, float r1, 
							float v10, float k1, float o1, float n1, 
							float d2, float S20, float r2, 
							float v20, float k2, float o2, float n2,
							float r, float *price, float *error) {


	// Indices used for trajectories and dimensions
	int ii; 
	// Indices needed to compute the sum end the sum square
	float sum, sum2;
	// The square root of the time increment
	float dt = sqrt((float)T/Nst);

	sum = 0.0;
	sum2 = 0.0;

	HestonSimulation_Modif (rnr, Ntraj, T, Nst,  
							S10, d1, v10, k1, o1, n1, r1,  
							S20, d2, v20, k2, o2, n2, r2, r);

    // payoff
	ComputePay(Str, Ntraj, Nst-1);

	for (ii = 0; ii < Ntraj; ii++) {
	   sum += TPE[ii];
	   sum2 += TPE[ii]*TPE[ii];
	}

	// Compute the price final error
	*price = exp(-rnr*T)*(sum/Ntraj);                               
 	*error = 1.96*sqrt((float)(exp(-2.0*rnr*T)/(Ntraj-1))*
			 (sum2 - (sum*sum)/Ntraj))/sqrt((float)Ntraj);

}



///////////////////////////////////////////////////////////////////////////////
// Some functions from http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c/
///////////////////////////////////////////////////////////////////////////////
static double MultModM (double a, double s, double c, double m)
   /* Compute (a*s + c) % m. m must be < 2^35.  Works also for s, c < 0 */
{
   double v;
   long a1;
   v = a * s + c;
   if ((v >= two53) || (v <= -two53)) {
      a1 = (long) (a / two17);
      a -= a1 * two17;
      v = a1 * s;
      a1 = (long) (v / m);
      v -= a1 * m;
      v = v * two17 + a * s + c;
   }
   a1 = (long) (v / m);
   if ((v -= a1 * m) < 0.0)
      return v += m;
   else
      return v;
}
/*-------------------------------------------------------------------------*/
static void MatVecModM (double A[3][3], double s[3], double v[3], double m)
   /* Returns v = A*s % m.  Assumes that -m < s[i] < m. */
   /* Works even if v = s. */
{
   int i;
   double x[3];
   for (i = 0; i < 3; ++i) {
      x[i] = MultModM (A[i][0], s[0], 0.0, m);
      x[i] = MultModM (A[i][1], s[1], x[i], m);
      x[i] = MultModM (A[i][2], s[2], x[i], m);
   }
   for (i = 0; i < 3; ++i)
      v[i] = x[i];
}
/*-------------------------------------------------------------------------*/
static void MatMatModM (double A[3][3], double B[3][3], double C[3][3],
                        double m)
   /* Returns C = A*B % m. Work even if A = C or B = C or A = B = C. */
{
   int i, j;
   double V[3], WA[3][3];
   for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j)
         V[j] = B[j][i];
      MatVecModM (A, V, V, m);
      for (j = 0; j < 3; ++j)
         WA[j][i] = V[j];
   }
   for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j)
         C[i][j] = WA[i][j];
   }
}
/*-------------------------------------------------------------------------*/
static void MatPowModM (double A[3][3], double B[3][3], double m, long n)
   /* Compute matrix B = A^n % m ;  works even if A = B */
{
   int i, j;
   double WA[3][3];

   /* initialize: W = A; B = I */
   for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; ++j) {
         WA[i][j] = A[i][j];
         B[i][j] = 0.0;
      }
   }
   for (j = 0; j < 3; ++j)
      B[j][j] = 1.0;

   /* Compute B = A^n % m using the binary decomposition of n */
   while (n > 0) {
      if (n % 2)
        MatMatModM (WA, B, B, m);
		MatMatModM (WA, WA, WA, m);
      n /= 2;
   }
}
/*-------------------------------------------------------------------------*/

/////////////////////////////////////////////////////////////////////////////
// Post initialization of the of the random number generator CMRG
/////////////////////////////////////////////////////////////////////////////
void PostInitCMRG(int Ntraj)
{
 int i;
 double s1[3] = {1, 1, 1};
 double s2[3] = {1, 1, 1};
 const int m1 = 2147483647;
 const int m2 = 2145483479;
 
 MatPowModM (A1, A1, m1, 1000000);
 MatPowModM (A2, A2, m2, 1000000);
 

 // Post-init seeds function of the process number
 for (i = 0; i < Ntraj; i++) {   	
	CMRG[0][i] = (int)s1[0];
	CMRG[1][i] = (int)s1[1];
	CMRG[2][i] = (int)s1[2];
	CMRG[3][i] = (int)s2[0];
	CMRG[4][i] = (int)s2[1];
	CMRG[5][i] = (int)s2[2];
	MatVecModM (A1, s1, s1, m1);
    MatVecModM (A2, s2, s2, m2);
 }	
}


static void dicho_fun_only (int Ntraj, int Nst, float rnr, float Str, float T, 
							float d1, float S10, float r1, float v10, float k1, 
							float o1, float n1, float d2, float S20, float r2, 
							float v20, float k2, float o2, float n2, float *pMin, 
							float *pMax, float *p_er, float *r, float *r_er, int NT){
	int i;
	float price, error;
	const float Mprice = *p_er;

	European_Heston(Ntraj, Nst, rnr, Str, T, d1, S10, r1, v10, k1, o1, n1, 
					d2, S20, r2, v20, k2, o2, n2, *r, &price, &error);

	for(i=1; i<NT; i++){
		if((Mprice > price)&&(Mprice < *pMax)){
			*r_er = *r_er/2.0f;
			*r -= *r_er;
			*pMin = price;
			European_Heston(Ntraj, Nst, rnr, Str, T, d1, S10, r1, v10, k1, o1, n1, 
							d2, S20, r2, v20, k2, o2, n2, *r, &price, &error); 
		}else{
			if((Mprice < price)&&(Mprice > *pMin)){
				*r_er = *r_er/2.0f;
				*r += *r_er;
				*pMax = price;
				European_Heston(Ntraj, Nst, rnr, Str, T, d1, S10, r1, v10, k1, o1, n1, 
								d2, S20, r2, v20, k2, o2, n2, *r, &price, &error);
			}
		} 
	}
	*p_er = error;

}

int main() {

float T;
float Str;
float rnr;
float S10, d1, v10, k1, o1, n1, r1;
float S20, d2, v20, k2, o2, n2, r2;
float M_price, pMax, pMin, p_er;
float r, r_er;
int YN;


int	Nst;						// Number of time steps
int Ntraj = 16384;				// Number of trajectories 

printf("Enter the maturity of the contract: ");
scanf("%f", &T);

// First test: on maturity
if (T>10.0){
	printf("This algorithm would require enormous time to converge \n");
	printf("Take a maturity smaller than 10 \n");
	exit(EXIT_FAILURE);
}

Nst = (int)(T/0.02);			// 50 time steps per unit of maturity		

printf("\n --> The payoff of the contract is given by: \n \n");
printf("(S1 - S2 + Str)_+, the strike 'Str' has to be positive \n");
printf("The contract should be at most 20%s ITM or OTM \n \n","%");
printf("Enter the positive value of S1_0: ");
scanf("%f", &S10);
printf("Enter the positive value of S2_0: ");
scanf("%f", &S20);
printf("Enter the value of the strike 'Str': ");
scanf("%f", &Str);

// Second test: on the payoff parameters
if((S10<=0)||(S20<=0)||(Str<0)){
	printf("Positive parameters are needed! \n");
	exit(EXIT_FAILURE);
}else{
	if(((S10+Str)/S20)>1.21){
		printf("The contract is too much ITM \n");
		exit(EXIT_FAILURE);
	}else{
		if((S20/(S10+Str))>1.21){
			printf("The contract is too much OTM \n");
			exit(EXIT_FAILURE);
		}
	}
}


printf("\n --> The model of each stock and its volatility is given by: \n \n");
printf("Stock 1:     dS1 = S1 (rnr - d1)dt + S1 sqrt(v1) dW1,  S1_0 = S10 \n");
printf("Volatility 1: dv1 = k1 (o1 - v1)dt + n1 sqrt(v1) dB1,  v1_0 = v10 \n \n");
printf("Stock 2:     dS2 = S2 (rnr - d2)dt + S2 sqrt(v2) dW2,  S2_0 = S20 \n");
printf("Volatility 2: dv2 = k2 (o2 - v2)dt + n2 sqrt(v2) dB2,  v2_0 = v20 \n \n");
printf("with: \n");
printf("<W1,B1>_t=r1*t, <W2,B2>_t=r2*t, <W1,W2>_t=r*t & <B1,B2>_t=r*r1*r2*t \n\n");


printf("Enter the value of n1:  ");
scanf("%f", &n1);
printf("Enter the value of k1:  ");
scanf("%f", &k1);
printf("Enter the value of o1:  ");
scanf("%f", &o1);
// Third test: on CIR parameters 
if((4*k1*o1)<(n1*n1)){
	printf("The relation: 4*k1*o1 >= n1*n1 must be satisfied \n \n");
	exit(EXIT_FAILURE);
}else{
// Fourth test: on CIR parameters considering the maturity 
	if(((2*k1*o1)<(n1*n1))&&(T>5)){
		printf("When the maturity is bigger than 5:\n");
		printf("the relation: 2*k1*o1 >= n1*n1 must be satisfied. \n");
		printf("Otherwise, the convergence is too slow! \n \n");
		exit(EXIT_FAILURE);
	}
}


printf("Enter the value of n2:  ");
scanf("%f", &n2);
printf("Enter the value of k2:  ");
scanf("%f", &k2);
printf("Enter the value of o2:  ");
scanf("%f", &o2);
// Third test: on CIR parameters 
if((4*k2*o2)<(n2*n2)){
	printf("The relation: 4*k2*o2 >= n2*n2 must be satisfied \n \n");
	exit(EXIT_FAILURE);
}else{
// Fourth test: on CIR parameters considering the maturity
	if(((2*k2*o2)<(n2*n2))&&(T>5)){
		printf("When the maturity is bigger than 5:\n");
		printf("the relation: 2*k2*o2 >= n2*n2 must be satisfied. \n");
		printf("Otherwise, the convergence is too slow! \n \n");
		exit(EXIT_FAILURE);
	}
}


printf("Enter the value of v10: ");
scanf("%f", &v10);
printf("Enter the value of v20: ");
scanf("%f", &v20);
// Fifth test: on the coherency between short-term and long-term volatility 
if((v10>3*o1)||(v20>3*o2)){
	printf("short-term volatility too big when compared to long-term one! \n");
	exit(EXIT_FAILURE);
}

printf("Enter the value of r1:  ");
scanf("%f", &r1);
printf("Enter the value of r2:  ");
scanf("%f", &r2);
// Sixth test: |r1| <= 1 and |r2| <= 1 
if((ABS(r1)>1)||(ABS(r2)>1)){
	printf("You gave either |r1| > 1 or |r2| > 1! \n");
	exit(EXIT_FAILURE);
}


printf("Enter the value of rnr: ");
scanf("%f", &rnr);
printf("Enter the value of d1:  ");
scanf("%f", &d1);
printf("Enter the value of d2:  ");
scanf("%f", &d2);

// Sixth test: |r1| <= 1 and |r2| <= 1 
if((ABS(r1)>1)||(ABS(r2)>1)){
	printf("You gave either |r1| > 1 or |r2| > 1! \n");
	exit(EXIT_FAILURE);
}


printf("Enter the market price of the contract: ");
scanf("%f", &M_price);

memory_allocation(Ntraj, Nst);
PostInitCMRG(Ntraj);


European_Heston(Ntraj, Nst, rnr, Str, T, d1, S10, r1, v10, k1, o1, n1, 
				d2, S20, r2, v20, k2, o2, n2, 0.9f, &pMin, &p_er);

European_Heston(Ntraj, Nst, rnr, Str, T, d1, S10, r1, v10, k1, o1, n1, 
				d2, S20, r2, v20, k2, o2, n2, -0.9f, &pMax, &p_er);

printf("\nThe minimum price, obtained for r=0.9, is: %f\n", (float)pMin);

printf("The maximum price, obtained for r=-0.9, is: %f\n", (float)pMax);

printf("The maximum error of the 95%% confidence interval = %f\n \n",(float)p_er);


// Seventh test: is the market price attainable ?
if((M_price>pMax)||(M_price<pMin)){
	if((M_price<(pMax+p_er))&&(M_price>pMax)){
		printf("The correlation r is in the interval [-1.0, -0.9] \n");
	}else{
		if((M_price>(pMin-p_er))&&(M_price<pMin)){
			printf("The correlation r is in the interval [0.9, 1.0] \n");
		}else{
			printf("The market price is not included in this interval, \n");
			printf("either there is no correlation for this market price or |r|~1. \n");
		}
	}
}else{
	printf("\nThis is a pure dichotomy algorithm, \n");
	r = 0.0f;
	r_er = 0.9f;
	p_er = M_price;
	dicho_fun_only (Ntraj, Nst, rnr, Str, T, d1, S10, r1, v10, k1, o1, n1, d2, 
					S20, r2, v20, k2, o2, n2, &pMin, &pMax, &p_er, &r, &r_er, 9);
	printf("Maximum price: %f \n",pMax);
	printf("Minimum price: %f \n",pMin);
	printf("The price error of the 95%% confidence interval = %f\n \n",(float)p_er);
	printf("Found correlation: %f \n",r);
	printf("Error on correlation: %f \n \n",r_er);
	printf("Is this precision sufficient? 0 for no ");
	scanf("%i", &YN);
	if(YN==0){
		printf("The number of trajectories is now equal to 65536 instead of 16384\n");
		Ntraj = 65536; 
		pMin -= p_er;
		pMax += p_er;
		p_er = M_price;
		free_memory(Nst);
		memory_allocation(Ntraj, Nst);
		PostInitCMRG(Ntraj);
		dicho_fun_only (Ntraj, Nst, rnr, Str, T, d1, S10, r1, v10, k1, o1, n1, d2, 
						S20, r2, v20, k2, o2, n2, &pMin, &pMax, &p_er, &r, &r_er, 5);
		printf("Maximum price: %f \n",pMax);
		printf("Minimum price: %f \n",pMin);
		printf("The price error of the 95%% confidence interval = %f\n \n",(float)p_er);
		printf("Found correlation: %f \n",r);
		printf("Error on correlation: %f \n \n",r_er);		
	}		
}

free_memory(Nst);
return 0;
}