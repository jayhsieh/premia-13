#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"pnl/pnl_vector.h"
#include"pnl/pnl_random.h"
#include"pnl/pnl_specfun.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_integration.h"

#define TAILLE_MAX 1000



typedef struct params {
	int n;
        double VarianceCurve;
        double Gamma;
        double Omega1;
        double Omega2;
        double Correl;
	

        double Small;
	double VIXFuture;
	double PutVIX;
	double Strike;
} params;

void Delete_params (params *p){	
	p->n=0;
	p->VarianceCurve=0;
	p->Gamma =0;
	 p->Omega1=0;
	p->Omega2=0;
	p->Correl=0;

	p->Small = 0;
	p->VIXFuture = 0;
	p->PutVIX = 0;
	p->Strike = 0;
}

 


double Hermit( double x, int n)
{
if (n<0)
return 0.0;
if (n ==0)
return x;
if(n==1)
return (x*x - 1.0);

return ( x * Hermit(x,n-1) - n * Hermit(x,n-2) );
}




double g(double x, params *p)
{
	int N ;
        double m_dVarianceCurve;
        double m_dGamma ;
        double m_dOmega1 ;
        double m_dOmega2 ;
        double m_dCorrel;
 
	N = p->n;
        m_dVarianceCurve= p->VarianceCurve;
        m_dGamma =  p->Gamma;
        m_dOmega1 = p->Omega1;
        m_dOmega2 = p->Omega2;
        m_dCorrel = p->Correl;

    return (m_dVarianceCurve * ((1.0 - m_dGamma) * exp(m_dOmega1 * x - m_dOmega1 * m_dOmega1 / 2.0) + m_dGamma * exp(m_dOmega2 * x - m_dOmega2 * m_dOmega2 / 2.0)));
 }

double z(double x, params *p)
{
	int N ;
        double m_dVarianceCurve;
        double m_dGamma ;
        double m_dOmega1 ;
        double m_dOmega2 ;
        double m_dCorrel;
 
	double LowerBond, UperBond;

	N = p->n;
        m_dVarianceCurve= p->VarianceCurve;
        m_dGamma =  p->Gamma;
        m_dOmega1 = p->Omega1;
        m_dOmega2 = p->Omega2;
        m_dCorrel = p->Correl;


                    if (m_dOmega1 == 0)
                            return 0;
                    
                    if (m_dOmega2 == 0)
                    {
                      
                            if (m_dVarianceCurve * m_dGamma == x)
                                return (m_dOmega1 / 2.0);
                            else
                                return (m_dOmega1 / 2.0 + log((x / m_dVarianceCurve - m_dGamma) / (1 - m_dGamma)) / m_dOmega1);
                       
                    }
                    if (x <= g(0,p))
                    {
                        UperBond = 0;
                        LowerBond = -1.0;
                        while (x <= g(LowerBond,p))
                        {
                            UperBond = LowerBond;
                            LowerBond -= 1.0;
                        }
                    }
                    else
                    {
                        UperBond = 1;
                        LowerBond = 0.0;
                        while (x >= g(UperBond,p))
                        {
                            LowerBond = UperBond;
                            UperBond += 1.0;
                        }
                    }
                    double error = 1, leght = 1;
                    double middle = (LowerBond + UperBond) / 2.0;
                    while (error >= p->Small && leght >= p->Small)
                    {
                        middle = (LowerBond + UperBond) / 2.0;
                        if (g(middle,p) < x)
                            LowerBond = middle;
                        else
                            UperBond = middle;
                        error = ABS(x - g(middle,p));
                        leght = UperBond - LowerBond;
                    }
                    return middle;

}

 double f_by_density(params *p)
        {
 
	int N ;
        double m_dVarianceCurve;
        double m_dGamma ;
        double m_dOmega1 ;
        double m_dOmega2 ;
        double m_dCorrel;

       int path,i ;
       double sup ;
       PnlVect *axis;

double res ;
	N = p->n;
        m_dVarianceCurve= p->VarianceCurve;
        m_dGamma =  p->Gamma;
        m_dOmega1 = p->Omega1;
        m_dOmega2 = p->Omega2;
        m_dCorrel = p->Correl;

       path = 2000;
       sup = 12;
       axis = pnl_vect_create(path);

            

            res = 0;
            
            for ( i = 0; i < path; i++)
                pnl_vect_set(axis,i, -sup + 2 * sup * i / (double)path);  
            
	    for ( i = 0; i < path - 1; i++)
	      res += (pnl_vect_get(axis,i + 1) - pnl_vect_get(axis,i)) *
                    sqrt(m_dVarianceCurve * (1 - m_dGamma) * exp(m_dOmega1 * pnl_vect_get(axis,i) - m_dOmega1 * m_dOmega1 / 2.0) +
                    m_dVarianceCurve * m_dGamma * exp(m_dOmega2 * pnl_vect_get(axis,i) - m_dOmega2 * m_dOmega2 / 2.0)) *
                    exp(-pnl_vect_get(axis,i) * pnl_vect_get(axis,i) / 2.0) / sqrt(2.0 *M_PI);
            return res;
        }
double call_by_density(double k, params *p)
        {
	int N ;
        double m_dVarianceCurve;
        double m_dGamma ;
        double m_dOmega1 ;
        double m_dOmega2 ;
        double m_dCorrel;

       int path,i ;
       double sup ;
       PnlVect *axis;
double strikeMin, res;

	N = p->n;
        m_dVarianceCurve= p->VarianceCurve;
        m_dGamma =  p->Gamma;
        m_dOmega1 = p->Omega1;
        m_dOmega2 = p->Omega2;
        m_dCorrel = p->Correl;

       path = 2000;
       sup = 12;
       axis = pnl_vect_create(path);

       strikeMin = z(k * k,p);
 

             res = 0;
            
            for (i = 0; i < path; i++)
	      pnl_vect_set(axis,i,strikeMin + sup * i / (double)path);
            for (i = 0; i < path - 1; i++)
	      res += (pnl_vect_get(axis,i+1) - pnl_vect_get(axis,i)) * MAX(sqrt(g(pnl_vect_get(axis,i),p)) - k, 0.0) *
                    exp(-pnl_vect_get(axis,i) * pnl_vect_get(axis,i) / 2.0) / sqrt(2.0 * M_PI);
            return res;
        }
double put_by_density(double k, params *p)
        {
	int N ;
        double m_dVarianceCurve;
        double m_dGamma ;
        double m_dOmega1 ;
        double m_dOmega2 ;
        double m_dCorrel;

       int path,i ;
       double inf ;
       PnlVect *axis;
	double strikeMax, res;

	N = p->n;
        m_dVarianceCurve= p->VarianceCurve;
        m_dGamma =  p->Gamma;
        m_dOmega1 = p->Omega1;
        m_dOmega2 = p->Omega2;
        m_dCorrel = p->Correl;

       path = 2000;
       inf = -12;
       axis = pnl_vect_create(path);

       strikeMax = z(k * k,p);
   
     
         res = 0;

            for (i = 0; i < path; i++)
	      pnl_vect_set(axis,i,inf + (i + 1) * (strikeMax - inf) / (double)path);
            for (i = 0; i < path - 1; i++)
	      res += (pnl_vect_get(axis,i+1) - pnl_vect_get(axis,i)) * MAX(k - sqrt(g(pnl_vect_get(axis,i),p)), 0.0) *
                    exp(-pnl_vect_get(axis,i) * pnl_vect_get(axis,i) / 2.0) / sqrt(2.0 * M_PI);
            return res;
        }


void Beta_Future_Put(double gamma, double zeta, params *p, double *betaF, double *betaP)
        {
            double put = 1;
            
	double temporary_gamma;
	double temporary_zeta;
	double temporary_beta;

double F, ErrorBeta, LeghtBeta, m_dVixPrice , m_dCallPrice, BetaMax , BetaMin;

            BetaMax = 1.0;
            BetaMin = 0.0;;


	temporary_gamma  =  gamma;
	temporary_zeta = zeta ;
	temporary_beta = 0.95;

            
          
	p->Gamma =temporary_gamma;
	 p->Omega1=temporary_zeta;
	p->Omega2=temporary_zeta*temporary_beta;
	
        F = f_by_density(p);
	ErrorBeta = 1.0;
        LeghtBeta = 1.0;


           p->Gamma = gamma;



          

	m_dVixPrice = p->VIXFuture  ;
	m_dCallPrice = p->PutVIX ;

           BetaMax = 1.0;
           BetaMin = 0.0;
            ErrorBeta = 1.0;
            LeghtBeta = 1.0;
            

            while (ErrorBeta >= 0.000001 && LeghtBeta > 0.000001)
            {
                LeghtBeta = -BetaMin + BetaMax;
                temporary_beta =(BetaMin + BetaMax) / 2.0;
                p->Omega2= temporary_zeta*temporary_beta;
                

                F = f_by_density(p);
                ErrorBeta = ABS(m_dVixPrice - F);
                if (F < m_dVixPrice)

                    BetaMax = temporary_beta;
                else
                    BetaMin = temporary_beta;



            }
             *betaF=temporary_beta;
            BetaMax = 1.0;
            BetaMin = 0.0;
            ErrorBeta = 1.0;
            LeghtBeta = 1.0;
            while (ErrorBeta >= 0.000001 && LeghtBeta > 0.000001)
            {
                LeghtBeta = -BetaMin + BetaMax;
                temporary_beta =(BetaMin + BetaMax) / 2.0;
                p->Omega2= temporary_zeta*temporary_beta;
               
                put = put_by_density(p->Strike,  p);

                if (put > m_dCallPrice)
                    BetaMax = temporary_beta;
                else
                    BetaMin = temporary_beta;

                ErrorBeta = ABS(put - m_dCallPrice);

            }
            *betaP =temporary_beta;
          


        }




void zetaMax(double gamma, params *p, double *omegaF, double *omegaP)
        {

            double put, F;
 


          
	double temporary_gamma;
	double temporary_zeta;
	double temporary_beta;

double omegaMinF, omegaMaxF, omegaMinP, omegaMaxP, ErrorOmega, LeghtOmega;

        ErrorOmega = 1.0;
        LeghtOmega = 1.0;

	temporary_gamma = gamma;
	temporary_zeta = 2.0;
	temporary_beta = 0.95;

	p->Gamma =temporary_gamma;
	 p->Omega1=temporary_zeta;
	p->Omega2=temporary_zeta*temporary_beta;

        omegaMinF = 0.0000000001; omegaMaxF = 15; omegaMinP = 0.000000000001; omegaMaxP = 12;

        ErrorOmega = 1.0;
        LeghtOmega = 1.0;

            while (ErrorOmega >= 0.00000001 && LeghtOmega > 0.000000001)
            {
                LeghtOmega = -omegaMinF + omegaMaxF;
                temporary_zeta = (omegaMinF + omegaMaxF) / 2.0;
                p->Omega1 =temporary_zeta;
                 p->Omega2 =0.0;
                
	F = f_by_density(p);

                ErrorOmega = ABS(p->VIXFuture  - F);
                if (F < p->VIXFuture )

                    omegaMaxF =temporary_zeta;
                else
                    omegaMinF =temporary_zeta;

            }
            *omegaF =temporary_zeta;
            ErrorOmega = 1.0;
            LeghtOmega = 1.0;
            while (ErrorOmega >= 0.00000001 && LeghtOmega > 0.000000001)
            {
                LeghtOmega = -omegaMinP + omegaMaxP;
                temporary_zeta = (omegaMinP + omegaMaxP) / 2.0;
                p->Omega1 =temporary_zeta;
                p->Omega2 =0.0;
                 put = put_by_density(p->Strike,  p);
                
                if (put >  p->PutVIX)
                    omegaMaxP = temporary_zeta;
                else
                    omegaMinP = temporary_zeta;

                ErrorOmega = ABS(put -  p->PutVIX);

            }
            *omegaP = temporary_zeta;
           


        }

void  gamma_ast(params *p, double *gamma)
        {
            double gammaMax = 1;
            double gamaMin = 0;
            double zetaF, zetaP;
            double error = 1;
            double lenght = 1;


            *gamma = (gammaMax + gamaMin) / 2.0; ;
            while (error > 0.0001 && lenght > 0.0001)
            {
                *gamma = (gammaMax + gamaMin) / 2.0;
                lenght = gammaMax - gamaMin;
                zetaMax((gammaMax + gamaMin) / 2.0 ,p, &zetaF, &zetaP);

                if (zetaF > zetaP)
                    gamaMin = (gammaMax + gamaMin) / 2.0;
                else
                    gammaMax = (gammaMax + gamaMin) / 2.0;
                error = ABS(zetaF - zetaP);

            }
	   

        }

void zeta_ast(double gamma,double *beta, double *zeta, params *p)
        {
            double zetaF_min, zetaP_min;
            double zetaF_max, zetaP_max;

double MinZeta , MaxZeta , error , lenght , betaF, betaP; 


	zetaMax(0.0,p,&zetaF_min, &zetaP_min);
	zetaMax(gamma,p,&zetaF_max, &zetaP_max);

            MinZeta =zetaF_min;
            MaxZeta =zetaF_max;
            error = 1;
            lenght = 1;
            
            *zeta = (MinZeta + MaxZeta) / 2.0;
            while (error > 0.000001 && lenght > 0.000001)
            {
                *zeta = (MinZeta + MaxZeta) / 2.0;
                lenght = MaxZeta - MinZeta;

		Beta_Future_Put(gamma, *zeta, p, &betaF, &betaP);

                
                if(betaF > betaP)
                    MinZeta = *zeta;
                else
                    MaxZeta = *zeta;
                error = ABS(betaF - betaP);

            }
*beta = betaF;
	   
           
        }

 void calibSkew(double k2, double Put2, params *p)
        {
	  int i;
            double gammaMin;
            double gammaMax;
            double aux, beta;
	 
int N ;
double  zetaF, zetaP,  Error, ErrorMethod ;
 
            int Index  ;
	gamma_ast(p, &gammaMin);
	gammaMax =MIN( p->VIXFuture * p->VIXFuture / p->VarianceCurve,p->Strike*p->Strike/p->VarianceCurve);

	    N = 10;
             
            
	PnlMat *point;
	PnlMat *skewmarket;
	point = pnl_mat_create( N,3);
	skewmarket = pnl_mat_create( N,2);
            
            
	pnl_mat_set(point,0,0,gammaMin);          


             zetaMax (gammaMin,p, &zetaF, &zetaP);
	pnl_mat_set(point,0,1,zetaF); 
	pnl_mat_set(point,0,2,0.0); 
            
           
            for ( i = 1; i < N; i++)
            {
                pnl_mat_set(point,i,0,gammaMin + i * (gammaMax - gammaMin) / (double)(N - 1));
                pnl_mat_set(skewmarket,i,0,pnl_mat_get(point,i,0));
 
zeta_ast(pnl_mat_get(point,i,0),&beta, &aux, p);
                pnl_mat_set(point,i,1,  aux  );
                
                 pnl_mat_set(point,i,2,  beta  ); 
                
            }

            Error = Put2, ErrorMethod = Put2;
            Index = 0;
            for( i = 0; i < N; i++)
            {
                p->Gamma = pnl_mat_get(point, i, 0);
                p->Omega1 = pnl_mat_get(point, i, 1);
                p->Omega2= pnl_mat_get(point, i, 1)*pnl_mat_get(point, i, 2);
                
		double putPrice =  put_by_density(k2,  p);
                ErrorMethod = ABS(putPrice - Put2);
                if (ErrorMethod < Error)
                {
                    Error = ErrorMethod;
                    Index = i;
                }
                
            }
            p->Gamma = pnl_mat_get(point, Index, 0);
            p->Omega1 = pnl_mat_get(point, Index, 1);
            p->Omega2= pnl_mat_get(point, Index, 1)*pnl_mat_get(point, Index, 2);
//printf("final Gamma = %f, Omega1 = %f, Omega2 = %f\n",p->Gamma, p->Omega1, p->Omega2);
        }


void RowFromFile(char *chaine, int numCol, PnlVect *res)
{
int i=0;
char delims[] = "\t";
char *result = NULL;
result = strtok( chaine, delims );
while( result != NULL ) {
pnl_vect_set(res,i, atof ( result ));

    result = strtok( NULL, delims );
i++;
}

}

double  calcul_var( double T,double Delta, double k1, double k2, double Theta , double RhoXY)
        {
            return ( SQR(1-Theta)*( exp( -2.0*k1*Delta) - exp( -2.0*k1*T) )/(2.0*k1) +  SQR(Theta)*(exp( -2.0*k2*Delta) - exp( -2.0*k2*T)  )/(2.0*k2)
                     + Theta*(1.0-Theta)*RhoXY* (exp( - (k1+k2)*Delta) - exp( -(k1+k2)*T)  )/(k1+k2)  );
          
             
        }


void calibToVIXSmiles(double k1, double k2, double Theta , double RhoXY, params *p, FILE* MarketData, FILE* OUTOMEGA)
{
char chaine[TAILLE_MAX] = "";
 int NumberMat,i,j;
PnlMat *Data;
PnlVect *aux;
PnlMat *OMEGA;
double varStateProces,T, Delta;

 if(MarketData != NULL)
    {
if(fgets(chaine, TAILLE_MAX, MarketData) != NULL)
i=0;
//LiborNumber = (int)atof ( chaine );
if(fgets(chaine, TAILLE_MAX, MarketData) != NULL)
NumberMat = (int) atof ( chaine );

Data = pnl_mat_create (NumberMat,7);
OMEGA = pnl_mat_create (2*NumberMat,3);
aux = pnl_vect_create (7);


for(j=0;j<NumberMat;j++)
{
	if(fgets(chaine, TAILLE_MAX, MarketData) != NULL) 
        	{
		 RowFromFile(chaine, NumberMat, aux);
		  for(i=0;i<7;i++)
		  { 
		   pnl_mat_set(Data,j,i,pnl_vect_get(aux,i)/100.0);
		   
		  }

        
		}
}

        fclose(MarketData);
    } 
 
for(i=0;i<NumberMat;i++)
{
p->VarianceCurve=pnl_mat_get(Data,i,1);
p->VIXFuture =pnl_mat_get(Data,i,2);
p->Strike  =pnl_mat_get(Data,i,3);
p->PutVIX  =pnl_mat_get(Data,i,4);

calibSkew(pnl_mat_get(Data,i,5), pnl_mat_get(Data,i,6), p);

if(i<NumberMat-1)
{
T = pnl_mat_get(Data,i+1,0);
Delta = pnl_mat_get(Data,i+1,0) - pnl_mat_get(Data,i,0);
}
else
{
Delta = pnl_mat_get(Data,i,0) - pnl_mat_get(Data,i-1,0);
T = pnl_mat_get(Data,i,0)+ Delta;

}
varStateProces = calcul_var(T, Delta,  k1,  k2,  Theta ,  RhoXY);
	
pnl_mat_set(OMEGA,2*i,0, T-Delta);
pnl_mat_set(OMEGA,2*i,1, p->Omega1);
 pnl_mat_set(OMEGA,2*i,2, p->VarianceCurve);


pnl_mat_set(OMEGA,2*i+1,0, T-p->Gamma*Delta);
pnl_mat_set(OMEGA,2*i+1,1, p->Omega2);
pnl_mat_set(OMEGA,2*i+1,2, p->VarianceCurve);

}
if (OUTOMEGA != NULL)
    {
fprintf( OUTOMEGA, "%d \n", 2*NumberMat ); // Ecriture du caractère A
for(i=0;i<2*NumberMat;i++)
        fprintf( OUTOMEGA, "%f \t %f \t %f \n", pnl_mat_get(OMEGA,i,0), pnl_mat_get(OMEGA,i,1), pnl_mat_get(OMEGA,i,2)); // Ecriture du caractère A
        fclose(OUTOMEGA);
    }

//pnl_mat_print(OMEGA);
}




int main()
{
 	
        double k1 ;
        double k2 ;
        double Theta ;
        double RhoXY ;
         
         
        
	params p;


	k1 = 0.35;
        k2 = 8.0;
        Theta = 0.35;
        RhoXY = 0.0;

	p.n = 0;
	p.Correl = 0;

	p.Small = 0.000001;
	 


 


FILE* Data = NULL;
Data = fopen("VIXDATA.txt", "r");

FILE* OUTOMEGA= NULL;
OUTOMEGA = fopen("OMEGA.txt", "w");


 


calibToVIXSmiles( k1,  k2,  Theta ,  RhoXY,&p, Data, OUTOMEGA);


return 0.0;
}

