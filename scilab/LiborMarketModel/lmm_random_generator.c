
#include"lmm_header.h"
#include "lmm_random_generator.h" 
#include <string.h>

#define lchaine 10000
#define INV_M (1.0/M)
#define INV_M1 (1.0/M1)

char ChBid[lchaine];
char *AL_ErrorMessage=ChBid;
char *Memory_Error="#Error : Memory allocation failure ";
char *Generator_Error="#Error : Unknown pseudo-random generator ";

const double pi=3.14159265358979;
const double DeuxPi=6.283185306;

static int FirstUse;


/* random generator algorithms */
static void KNUTH(double *sample);
static void MRGK3(double *sample);
static void MRGK5(double *sample);
static void SHUFL(double *sample);
static void LECUYER(double *sample);
static void TAUS(double *sample);

static void (*CurrentGenerator)(double *sample);
static void InitGenerator(char *ErrorMessage, char *Generator_Name);

double Generator();
int Uniform(int max);
double Gauss();


int mallocRandomGenerator(RandomGenerator **ptRand ,int numOfFac )
{
  int i;
  RandomGenerator *pt;
  
  InitGenerator(AL_ErrorMessage,"LECUYER");
  if (strcmp("",AL_ErrorMessage)!=0) 
    {
      printf("%s\n",AL_ErrorMessage);
      return;
    }
  
  pt=(RandomGenerator *)malloc(sizeof(RandomGenerator));

  pt->numberOfFactors=numOfFac;
  pt->val=(double*)malloc(sizeof(double)*pt->numberOfFactors);
  randomVector(pt);
  *ptRand=pt;
  return(1);
}

int freeRandomGenerator(RandomGenerator **ptRand)
{
  int i;
  RandomGenerator  *pt;

  pt=(*ptRand);
  *ptRand=NULL;
  free(pt->val);
  return(1);
}

int randomVector(RandomGenerator *ptRand)
{// compute a new random gaussian vector 
  int i;
  
  for (i=0; i<ptRand->numberOfFactors; i++)
    {
      ptRand->val[i]=Gauss();
    }
};


int printRandomGenerator(RandomGenerator *ptRand)
{
  int i;
  for(i=0;i<ptRand->numberOfFactors;i++)
    {
      printf("%lf \n",ptRand->val[i]);
      
    }
  printf("\n");
}
  

double getRandom(RandomGenerator *ptRand, int factorNumber)
{
  if (factorNumber> ptRand->numberOfFactors)
    {
      printf("factor number out of range \n");
      exit(1);
    }
  else
    {
      return(ptRand->val[factorNumber]);
    }
  
}


/*********************************************************************************************************************************************************
 *
 *
 *
 *   Interface with random generators
 *
 *
 *
 *
 *
 *
 ******************************************************************************************************************************************************/

double Generator()
{
  static double aux;
  CurrentGenerator(&aux);
  return aux;
}

int Uniform(int max)
{
  double aux;
  /*uniform law on {1,...,max}*/
  CurrentGenerator(&aux);
  return (int)(aux*(double)max);
}

double Gauss()
{
  double aux1,aux2;
  /*BoxMuller gaussian generator*/
  CurrentGenerator(&aux1);
  CurrentGenerator(&aux2);
  return sqrt(-2*log(aux1))*cos(DeuxPi*aux2);
}
  


/*****************************************************************************************************************************************************
 *
 *
 *  Random generator algorithms from earlier version of premia
 *
 *
 ****************************************************************************************************************************************************/

static void KNUTH(double *sample)
{
  static long M= 1000000000;
  static long SEED= 161803398;

  /* Initialize the sequence with a positive seed */
  static int alea= 1;
  static int inc1, inc2;
  static long t_alea[56];
  long X_n, y_k;
  int i, ii, l;
  
  /* First call to the sequence */
  if (FirstUse)
    {
      X_n= SEED- alea;
      X_n%= M;
      t_alea[55]= X_n;
      y_k= 1;
      /* Initialization of the table */
      for(i= 1; i<= 54; i++)
		{
		  ii= (21*i)%55; /* 21 was chosen to alleviate initial
							nonrandomness problems */
		  t_alea[ii]= y_k;
		  y_k= X_n - y_k;
		  if(y_k < 0)
			y_k+= M;
		  X_n= t_alea[ii];
		}
      
      /* Randomization of the elements of the table */
      for(l=  1; l<= 4; l++)
		{
		  for(i= 1; i<= 55; i++)
			{
			  t_alea[i]-= t_alea[1+(i+30)%55];
			  if(t_alea[i] < 0)
				t_alea[i]+= M;
			}
		}
      inc1= 0;
      inc2= 31;  /* 31 is a special value of Knuth : 31= 55-24 */
      alea= 1;
	  FirstUse=0;
    }
  
  /* For each call to the sequence, computation of a new point */
  if(++inc1 == 56)
    inc1= 1;
  if(++inc2 == 56)
    inc2= 1;
  /* Substractive method*/
  X_n= t_alea[inc1] - t_alea[inc2];
  
  if(X_n < 0)
    X_n+= M;
  t_alea[inc1]= X_n;
  /* Normalized value */
  *sample=X_n*INV_M;
  return;
}

/* ----------------------------------------------------------------------- */
/* Combination of two multiplicative recursive generators of order 3 (k=3) */
/* ----------------------------------------------------------------------- */
static void MRGK3(double *sample)
{
  static double M1= 4294967087.0;
  static double M2= 4294944443.0;
  static double A12= 1403580.0;
  static double A13N= 810728.0;
  static double A21= 527612.0;
  static double A23N= 1370589.0;
  static double NORM= 2.328306549295728e-10;

  static double x10, x11, x12, x20, x21, x22;
  long k;
  double p1, p2;
  
  /* First call to the sequence */
  if (FirstUse)
    {
      /* Initialization */
      x10=231458761.;
      x11=34125679.;
      x12=45678213.;
      x20=57964412.;
      x21=12365487.;
      x22=77221456.;
      FirstUse=0;
    }
  
  /* For each call to the sequence, computation of a new point */
  /* First generator */
  p1= A12*x11 - A13N*x10;
  k= (long)floor(p1/M1); /*TOCHECK*/
  p1-= k*M1;

  if(p1 < 0.0)
    p1+= M1;

  x10= x11;
  x11= x12;
  x12= p1;
  
  /* Second generator */
  p2= A21*x22 - A23N*x20;
  k= (long)floor(p2/M2);/*TOCHECK*/
  p2-= k*M2;

  if(p2 < 0.0)
    p2+= M2;

  x20= x21;
  x21= x22;
  x22= p2;

  /* Combination of the two generators */
  if (p1< p2)
    *sample= (p1- p2+ M1)*NORM;
  else
    *sample=(p1- p2)*NORM;
  return;
}

/* ----------------------------------------------------------------------- */
/* Combination of two multiplicative recursive generators of order 5 (k=5) */
/* ----------------------------------------------------------------------- */
static void MRGK5(double *sample)
{
  static double M1= 4294949027.0;
  static double M2=    4294934327.0;
  static double A12=   1154721.0;
  static double A14=   1739991.0;
  static double A15N=  1108499.0;
  static double A21=   1776413.0;
  static double A23=   865203.0;
  static double A25N= 1641052.0;
  static double NORM=  2.3283163396834613e-10;
  
  static double x10, x11, x12, x13, x14, x20, x21, x22, x23, x24;
  long k;
  double p1, p2;
  
  /* First call to the sequence */
  if (FirstUse)
    {
      /*Initialization*/
      x10= 231458761.;
      x11= 34125679.;
      x12= 45678213.;
      x13= 7438902.;
      x14= 957345.;

      x20= 57964412.;
      x21= 12365487.;
      x22= 77221456.;
      x23= 816403.;
      x24= 8488912.;
      FirstUse=0;
    }
  
  /* For each call to the sequence, computation of a new point */
  /* First generator with Schrage method */
  p1= A12*x13 - A15N*x10;
  
  if(p1> 0.0)
    p1-= A14*M1;

  p1+= A14*x11;
  k= (long)floor(p1/M1);/*TOCHECK*/
  p1-= k*M1;

  if(p1< 0.0)
    p1+= M1;
  
  x10= x11;
  x11= x12;
  x12= x13;
  x13= x14;
  x14= p1;
  
  /* Second generator with Schrage method */
  p2= A21*x24 - A25N*x20;

  if(p2> 0.0)
    p2-= A23*M2;

  p2+= A23*x22;
  k= (long)floor(p2/M2);/*TOCHECK*/
  p2-= k*M2;

  if(p2< 0.0)
    p2+= M2;

  x20= x21;
  x21= x22;
  x22= x23;
  x23= x24;
  x24= p2;

  /*Combination of the two generators */
  if (p1<= p2)
    *sample= (p1- p2+ M1)*NORM;
  else
    *sample=(p1- p2)*NORM;

  return;
}
/* ------------------------------------------------------------------------ */
/* Random numbers generator of  Park & Miller with Bayes & Durham shuffling
   procedure : the next random number is not obtained from the previous one
   but we use an intermediate table which contains the 32 precedent random
   numbers and we choose one of them randomly.  */
/* ----------------------------------------------------------------------- */
static void SHUFL(double *sample)
{
  static long A= 16807;        /* multiplier */
  static long M= 2147483647L;    /* 2**31 - 1  */
  static long Q= 127773L;        /*   M div A  */
  static long R= 2836;           /*   M mod A  */
  long N1;
  int j;
  long hi;                 /* high order bit */
  static long y= 0;
  static long t[32];       /* 32 refers to the size of a computer word */
  /* Initialisation */
  static long x;
  
  N1=(M/32);
  
  /* First call to the sequence */
  if (FirstUse)
    {
      x= 1043618065;
      
      /* After 8 "warm-ups", initialisation of the shuffle table */
      for (j= 39; j>= 0; j--)
		{
		  hi= x/Q;
		  /*Schrage's method to avoid overflows */
		  x= A*(x- hi*Q)- R*hi;
		  if (x < 0)
			x+= M;
		  if (j< 32)
			t[j]= x;
		}
      y= t[0];
	  FirstUse=0;
    }
  
  /* For each call to the sequence, computation of a new point */
  hi= x/Q;
  x= A*(x-hi*Q)- R*hi;
  if (x < 0)
    x+= M;

  /* Shuffling procedure of Bayes & Durham */
  /* Index j dependent on the last point */
  j= y/N1;
  /* Next point dependent on j */
  y= t[j];
  
  t[j]= x;
  *sample= INV_M*y;
  
  return;
}

/* ------------------------------------------------------------------------- */
/* Random numbers generator of L'Ecuyer with Bayes & Durham shuffling
   procedure :
   Combination of two short periods LCG to obtain a longer period generator.
   The period is the least common multiple of the 2 others.  */
/* ------------------------------------------------------------------------- */

static void LECUYER(double *sample)
{
  static long A1= 40014;        /* multiplier of the 1st generator */
  static long A2= 40692;        /* multiplier of the 2nd generator */
  static long M1= 2147483647;   /* 2**31 - 1   */
  static long M2= 2147483399;   /* 2**31 - 249 */
  static long Q1= 53668;        /* m1 div a1   */
  static long Q2= 52774;        /* m2 div a2   */
  static long R1= 12221;        /* m1 mod a1   */
  static long R2= 3791;         /* m2 mod a2   */
  long N1;

  static long x;
  static long y= 978543162;
  int j;
  long hi;               /* high order bit */
  static long z= 0;
  static long t[32];     /* 32 is the size of a computer word */

  N1= (M1/32);
  
  /* First call to the sequence */
  if (FirstUse)
    {
      x= 437651926;
      y= x;
      
      /* After 8 "warm-ups", initialisation of the shuffle table */
      for (j= 39; j>= 0; j--)
		{
		  /* Park & Miller's generator */
		  hi= x/Q1;
		  x= A1*(x-hi*Q1) - R1*hi;
		  if (x < 0)
			x+= M1;
		  if (j < 32)
			t[j]= x;
		}
      z= t[0];
	  FirstUse=0;
    }
  
  /* For each call to the sequence, computation of a new point */
  /* First generator */
  hi= x/Q1;
  x= A1*(x-hi*Q1) - R1*hi;
  if (x < 0)
    x+= M1;

  /* Second generator */
  hi= y/Q2;
  y= A2*(y-hi*Q2) - R2*hi;
  if (y < 0)
    y+= M2;
  
  /* Shuffling procedure of Bayes & Durham */
  /* Index j dependent on the last point */
  j= z/N1;
  /* Next point dependent on j */
  z= t[j]- y;
  t[j]= x;
  
  
  /* To avoid 0 value */
  if (z < 1)
    z+= M1-1;
  
  *sample= INV_M1*z;
  
  return;
}

/* ------------------------ */
/* Tausworthe Algorithm   */
/* ------------------------ */

/* maximal number of combined generators */
#define TAUS_MAX 10

/* ---------------------------------------------------- */
/* Generation of a random bit
   Algorithm based on a prime polynomial :
   here we choose x^18 + x^5 + x^2 + x + 1 .
   It is described in 'Numerical Recipes in C' page 296. */
/* ---------------------------------------------------- */
int bit_random()
{
  static int compt = 1;
  int degre = 18;
  static unsigned long a;
  unsigned long new_bit;

  /* Initialisation for the n first values */
  /* random number over [1, 2^18] ; 2^18= 262144 */
  if(compt == 1)
    {
      a= 176355;
    }
  compt++;
    
  /* Next bit calculation by the recurrence relation  */
  new_bit= (a & (1<<17)) >> 17
    ^ (a & (1<<4)) >> 4
    ^ (a & (1<<1)) >> 1
    ^ (a & 1);
  a <<= 1;
  /* The new bit is shift on the right */
  a ^= new_bit;
  /* The most left bit is put to 0 */
  a ^= (1 << degre);
  
  return((int)new_bit);
}

/* --------------------------------------------------------------- */
/* Generation of a word of k random bits. */
/* --------------------------------------------------------------- */
unsigned long random_word(int k)
{
  int i, bit;
  unsigned long mot;
  
  mot= 0;
  for(i= 0; i< k; i++)
    {
      bit= bit_random();
      mot= (mot<<1) ^ bit;
    }
  return(mot);
}

/* ---------------------------------------------------------------- */
/*  Tausworthe  Algorithm
   Combination  of J Tausworthe generators
             u(n)[j]= u(n-r)[j] ^ u(n-k)[j],
   with parameters k, q, r, s and t.
   Generator :
             v= =(u[0] ^ u[1] ... ^ u{J-1])/2^32.
     
   L= 32 length of a word.    */
/* ---------------------------------------------------------------- */
static void TAUS(double *sample)
{
  int L= 32;
  int J= 3;
  
  int i;
  static unsigned long u[TAUS_MAX];
  static unsigned long c[TAUS_MAX];
  unsigned long b;
  static int k[TAUS_MAX], q[TAUS_MAX];
  static int s[TAUS_MAX], r[TAUS_MAX], t[TAUS_MAX];
  unsigned long v= 0;

  /* First call to the sequence. Initialisation */
  if (FirstUse)
    {
      /* Choice of the parameters to have ME-CF generators (cf L'Ecuyer) */
      k[0]= 31; q[0]= 13; s[0]= 12;
      r[0]= k[0]- q[0]; t[0]= k[0]- s[0];
      
      k[1]= 29; q[1]= 2; s[1]= 4;
      r[1]= k[1]- q[1]; t[1]= k[1]- s[1];
      
      k[2]= 28; q[2]= 3; s[2]= 17;
      r[2]= k[2]- q[2]; t[2]= k[2]-s[2];
      
      /* constant c : k bits to one and (L-k) bits to zero */
      /* c[j]= 2^32 - 2^(L-k[j]) */
      c[0]= 4294967294ul;
      c[1]= 4294967288ul;
      c[2]= 4294967280ul;
      
      /* initialisation of each generator */
      u[0]= 0; u[1]= 0; u[2]= 0;
      for(i= 0; i< J; i++)
		{
		  /* The first k bits are chosen randomly */
		  u[i]= random_word(k[i]);
		  /* The next L-k bits are initially fixed to zero */
		  u[i] <<= (L- k[i]);
		  /* Now they are computed with the recurrence on u */
		  b= u[i] << q[i];
		  b ^= u[i];
		  b >>= k[i];
		  u[i] ^= b;
		}
      FirstUse=0;
    }
  
  /* For each call to the sequence, computation of a new point */
  for(i= 0; i< J; i++)
    {
      /* Calculus of the next point for the J generators */
      /* 6 steps explained by L'Ecuyer */
      b=(u[i]<<q[i]) ^ u[i];         /* Steps 1 and 2 */
      b >>= t[i];                    /* Step 3 */
      u[i]= (u[i] & c[i]) << s[i];   /* Steps 4 et 5 */
      u[i] ^= b;                     /* Step 6 */
      /* Combination : XOR between the J generators */
      v^= u[i];
    }

  /* Normalization by 1/2^32 */
  *sample=v * 2.3283064365e-10;

  return;
}

static void InitGenerator(char *ErrorMessage, char *Generator_Name)
{
  /*initialization of the current generator*/
  if (strcmp("KNUTH",Generator_Name)==0) CurrentGenerator=KNUTH;
  else if (strcmp("MRGK3",Generator_Name)==0) CurrentGenerator=MRGK3;
  else if (strcmp("MRGK5",Generator_Name)==0) CurrentGenerator=MRGK5;
  else if (strcmp("SHUFL",Generator_Name)==0) CurrentGenerator=SHUFL;
  else if (strcmp("LECUYER",Generator_Name)==0) CurrentGenerator=LECUYER;
  else if (strcmp("TAUS",Generator_Name)==0) CurrentGenerator=TAUS;
  else {
	strcat(ErrorMessage,Generator_Error);
	strcat(ErrorMessage,"(");
	strcat(ErrorMessage,Generator_Name);
	strcat(ErrorMessage,")");
  }
  FirstUse=1;
}


