/*************************************************************************************************************************************
 *
 *
 *
 *   Simulation in the libor market model of a  swaption
 *          method:  Bushy Tree
 *
 *   Marouen Messaoud
 *   March  2004
 *
 *
 *
 *
 ************************************************************************************************************************************/

#include "MMsimulBGM.h"
#include "spline.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "optim-code/QuasiNewton.h"
void matriceprod(double ** A,double **B,double **res,long nl, long nc)
{

	int k,i,j,m;
	double **temp;
	temp = (double **)malloc((nl)*sizeof(double*));
	for(m=0;m<nl;m++)
	{
		temp[m] = (double *)malloc(nc*sizeof(double));
		for(j=0;j<nc;j++)
		{
			temp[m][j] =0.0;
		}
	}

	for(i=0;i<nl;i++)
	{
		for(j=0;j<nc;j++)
		{
			for(k=0;k<nc;k++)
			{
				temp[i][j]+= A[i][k]*B[k][j];
			}
		}
	}

	for(i=0;i<nl;i++)
	{
		for(j=0;j<nc;j++)
		{
			res[i][j]=temp[i][j];
		}
	}
	free(temp);
}

double costFunction ( double *lamda )
{
	double res,***R,**Stilde,**STemp;
	int k,l,i,j,size,m,z,h;
	int NFactor;
	NFactor = NBranches-1;
	res = 0.0;
	size = (NFactor)*(NFactor-1)/2;
	
	Stilde = (double **)malloc(NFactor*sizeof(double*));
	for(m=0;m<NFactor;m++)
	{
		Stilde[m] = (double *)malloc(NFactor*sizeof(double));
	}

	R = (double ***) malloc(size*sizeof(double**));
	for(m=0;m<size;m++)
	{
		R[m] = (double **) malloc(NFactor*sizeof(double**));
		for(i=0;i<NFactor;i++)
		{
			R[m][i] = (double *) malloc(NFactor*sizeof(double));
			for(j=0;j<NFactor;j++)
			{
				R[m][i][j] =0.0;
			}
		}
	}
	m=0;
	for(k=1;k<NFactor;k++)
	{
		for(l=k+1;l<NFactor+1;l++)
		{
			for(i=0;i<NFactor;i++)
			{
				R[m][i][i]=1;
			}
			R[m][k-1][l-1] = pow(-1,m)*sin(lamda[m]);
			R[m][l-1][k-1] = pow(-1,m+1)*sin(lamda[m]);
			R[m][k-1][k-1] = cos(lamda[m]);
			R[m][l-1][l-1] = cos(lamda[m]);
			m++;			
		}
	}
	
	for(m=0;m<NFactor;m++)
	{
		for(j=0;j<NFactor;j++)
		{
			Stilde[m][j] =R[0][m][j];
		}
	}
	for(k=0;k<size-1;k++)
	{
		matriceprod(Stilde,R[k+1],Stilde,NFactor, NFactor);
	}

	STemp=(double **) malloc(NBranches*sizeof(double));
	for (h=0;h<NBranches;h++){
		STemp[h] =(double*) malloc((NBranches-1)*sizeof(double));
	}
	matriceprod(S,Stilde,STemp,NFactor+1, NFactor);
	res = 0.0;
	z= (NFactor+1)/2;
	for(i=0;i<(NFactor+1)/2;i++)
	{
		for(j=0;j<z;j++)
		{
			res += pow(STemp[i][j]+STemp[NFactor-i][j],2);
		}
	}
	if(flag==1)
	{
		for(i=0;i<(NBranches);i++)
		{
			for(j=0;j<NBranches-1;j++)
			{
				Mattemp[i][j] = STemp[i][j];
			}
		}
	}
	if (flag==2)
	{
		for(i=0;i<(NBranches);i++)
		{
			for(j=0;j<NBranches-1;j++)
			{
				S[i][j] = STemp[i][j];
			}
		}
	}

	for(m=0;m<size;m++)
	{
		for(i=0;i<NFactor;i++)
		{
			free(R[m][i]);

		}
		free(R[m]);
	}
	free(R);

	for(m=0;m<NFactor;m++)
	{
		free(Stilde[m]);
	}
	free(Stilde);
	for (h=0;h<NBranches;h++){
		free(STemp[h]);
		STemp[h]=NULL;
	}
	free(STemp);
	STemp=NULL;
	return res;
}

void gradCostFunction ( double *lamda , double *grad )
{

	int p;
	double res,***R,**Stilde,cost,**STemp;
	int k,l,i,j,size,m,z,h;
	int NFactor;
	NFactor = NBranches-1;
	size = (NFactor)*(NFactor-1)/2;
	//z= (size+1)/2;
	flag=1;
	cost = costFunction (lamda );
	flag=0;

	STemp=(double **) malloc(NBranches*sizeof(double));
	for (h=0;h<NBranches;h++){
		STemp[h] =(double*) malloc((NBranches-1)*sizeof(double));
	}

	Stilde = (double **)malloc(NFactor*sizeof(double*));
	for(m=0;m<NFactor;m++)
	{
		Stilde[m] = (double *)malloc(NFactor*sizeof(double));
	}

	R = (double ***) malloc(size*sizeof(double**));
	for(m=0;m<size;m++)
	{
		R[m] = (double **) malloc(NFactor*sizeof(double**));
		for(i=0;i<NFactor;i++)
		{
			R[m][i] = (double *) malloc(NFactor*sizeof(double));
			for(j=0;j<NFactor;j++)
			{
				R[m][i][j] =0.0;
			}
		}
	}

	for(p=0;p<size;p++)
	{

		m=0;
		for(k=1;k<NFactor;k++)
		{
			for(l=k+1;l<NFactor+1;l++)
			{
				if(m==p)
				{
					for(i=0;i<NFactor;i++)
					{
						R[m][i][i]=0.0;
					}
					R[m][k-1][l-1] = pow(-1,m)*cos(lamda[m]);
					R[m][l-1][k-1] = pow(-1,m+1)*cos(lamda[m]);
					R[m][k-1][k-1] = -sin(lamda[m]);
					R[m][l-1][l-1] = -sin(lamda[m]);
					m++;
				}else{
					for(i=0;i<NFactor;i++)
					{
						R[m][i][i]=1;
					}
					R[m][k-1][l-1] = pow(-1,m)*sin(lamda[m]);
					R[m][l-1][k-1] = pow(-1,m+1)*sin(lamda[m]);
					R[m][k-1][k-1] = cos(lamda[m]);
					R[m][l-1][l-1] = cos(lamda[m]);
					m++;
				}
			}
		}

		for(m=0;m<NFactor;m++)
		{
			for(j=0;j<NFactor;j++)
			{
				Stilde[m][j] = R[0][m][j];
			}
		}
		for(k=0;k<size-1;k++)
		{
			matriceprod(Stilde,R[k+1],Stilde,NFactor, NFactor);
		}
		matriceprod(S,Stilde,STemp,NFactor+1, NFactor);
		res = 0.0;
		z= (NFactor+1)/2;
		for(i=0;i<(NFactor+1)/2;i++)
		{
			for(j=0;j<z;j++)
			{
				res += 2*(Mattemp[i][j]+Mattemp[NFactor-i][j])*(STemp[i][j]+STemp[NFactor-i][j]) ;
			}
		}

		grad[p]= res;
	}
	for(m=0;m<size;m++)
	{
		for(i=0;i<NFactor;i++)
		{
			free(R[m][i]);
			
		}
		free(R[m]);
	}
	free(R);
	for(m=0;m<NFactor;m++)
	{
		free(Stilde[m]);
	}
	free(Stilde);
	for (h=0;h<NBranches;h++){
		free(STemp[h]);
		STemp[h]=NULL;
	}
	free(STemp);
	STemp=NULL;
}
int RotateSimplex(int NFactor)
{
	double *lamda,*grad,*sol;
	double grad_tol , step_tol ,lambda0;
	long size;
	int h,max_counter,verbosity,saveSuccessiveXinFile;

	Mattemp = (double **) malloc(NBranches*sizeof(double));
	for (h=0;h<NBranches;h++){
		Mattemp[h] =(double*) malloc((NBranches-1)*sizeof(double));
	}
	flag=0;
	max_counter=100;
	verbosity=0;
	saveSuccessiveXinFile=0;
	grad_tol = 0.00000000000001;
	step_tol = 0.00000000000001;
	lambda0 = 1.325447;
	size = (NFactor*(NFactor-1)/2);
	lamda =(double *)malloc(size*sizeof(double));

	grad =(double *)malloc(size*sizeof(double)) ;
	sol =(double *)malloc(size*sizeof(double)) ;
	QuasiNewton( size ,
	                lamda,
	                costFunction,
	               	gradCostFunction,
	                sol,
	                grad_tol,
	                step_tol,
	                max_counter,
	                verbosity,
	                saveSuccessiveXinFile);
	flag=2;
	costFunction(sol);
	flag=0;
	for (h=0;h<NBranches;h++){
		free(Mattemp[h]);
		Mattemp[h]=NULL;
	}

	free(Mattemp);
	Mattemp=NULL;
	free(lamda);
	free(sol);
	free(grad);
	return 0;
}
int SetS(){
//S a une taille NbFactor=NBranches-1 colonnes et NBranches lignes
	long kl,kc;
	kl=0;
	for(kc=kl;kc<NBranches-1;kc++){
		S[kl][kc]=-sqrt( NBranches/((kc+2.0)*(kc+1.0)) );
	}
	for(kl=1;kl<NBranches;kl++){
		for(kc=0;kc<kl-1;kc++){
			S[kl][kc]=0.0;
		}
		kc=kl-1;
		S[kl][kc]=sqrt( NBranches*(kc+1.0)/((kc+2.0)) );
		for(kc=kl;kc<NBranches-1;kc++){
			S[kl][kc]=-sqrt( NBranches/((kc+2.0)*(kc+1.0)) );
		}
	}

	return 0;
}
int setCovarMatrix(double **C ,long h)
{
	// il faut construire la matrice de covariance a la date h
	double a, b, c, d, beta,Ti,Tj;
	int i,j;
	a=-0.02;
	b=0.5;
	c=1;
	d= 0.1;
	beta= 0.1;
	for (i=0;i<NRates;i++){
		Ti=Tau[NSteps] +i*period;
		for(j=0;j<NRates;j++)
		{
			Tj=Tau[NSteps] +j*period;
			C[i][j]=(( a+b*(Ti-Tau[h]) )*exp(-c*(Ti-Tau[h])) + d)* ( ( a+b*(Tj-Tau[h]) )*exp(-c*(Tj-Tau[h]))+d)*exp(-beta*fabs(Tj-Ti));// covariance croissante en temps!!
 			
		}
	}
	return 0;
}

int SetBbar(long h)
{
	long i,j,k;
	for(i=0;i<NRates;i++)
	{
		for(j=0;j<NBranches;j++)
		{
			for(k=0;k<NBranches-1;k++)
			{
				Bbar[h][i][j]+=A[i][k]*S[j][k];
			}
		}
	}
	return 0;
}
int SetLogShiftOfBranch(long h){
	long i,k;
	double tauh,sqrttauh,Cii;
	tauh=Tau[h+1]-Tau[h];
	sqrttauh=sqrt(tauh);
	for (i=0;i<NRates;i++){
		Cii = C[h][i][i];
		for(k=0;k<NBranches;k++){
			LogShiftOfBranch[h][k][i] =-0.5*Cii*tauh+sqrttauh*Bbar[h][i][k];

		}
	}

	return 0;
}

/* return the payoff at time h*/
double Intrinsic(long h){
	// to fill and get the terminal intrinsec value on the tree
	return 0.0;
}
double CheckForEarlyExercise(long h,double esp )
{
	double res=esp*h;
	// to fill if there is a earlyexercice on the tree
	return res;
}
//M as Input : the covariance matrix at time t to t+dt
//C as output : the decorrelated covariance matrix C(i,i)=M(i,i)
int DecorrelateCovarianceMatrix(double **M, int NRate,int NFactor, double ** C,int h)
{
	int i,j,k;
	double somme=0.0,sqrtTemp;
	double *Mdiag=(double *) malloc(NRate*sizeof(double));

	for( i=0;i<NRate;i++)
	{
		Mdiag[i]=M[i][i];
	}
	//AfficheMat(M,NRate,NRate);
	Cholesky(M, NRate );

	for(i=0;i<NRate;i++)
	{
		somme=0.0;
		for(k=0; k<NFactor && k<=i ;k++)
		{
			somme+=M[i][k]*M[i][k];
		}
		sqrtTemp= sqrt(Mdiag[i]/somme);
		for(j=0;j<=i && j<NFactor ;j++)
		{
			A[i][j]=M[i][j]*sqrtTemp;
		}
		for(j=i+1;j<NFactor;j++)
		{
			A[i][j]=0.0;
		}
	}
	//AfficheMat(A,NRate,NFactor);
	//C=A*A';

	for(i=0;i<NRate;i++)
	{
		for(j=0;j<NRate;j++)
		{
			C[i][j]=0.0;
		}
	}

	for(j=0;j<NRate;j++)
	{
		for( k=0;k<NFactor;k++)
		{
			if(k==0)
				//Vol.volElt(j,k,i) = 0.2;
				A[j][k] = 0.2;
			else
			{
				//Vol.volElt(j,k,i) = 0.1 + 0.07*((double)i/(number_of_period-1));
				A[j][k] = 0.1 + 0.07*((double)h/(NSteps-1));
			}
			//printf(" vol i %d j %d k %d val %f \n", h , j , k ,A[j][k] );
		}

	}


	for(i=0;i<NRate;i++)
	{
		for(j=0;j<=i;j++)
		{
			for(k=0 ; k<NFactor;k++)
			{
				C[i][j]+=A[i][k]*A[j][k];
			}
			C[j][i]=C[i][j];
		}

	}
	//AfficheMat(C,NRate,NRate);
	free( Mdiag);
	Mdiag =NULL;

	return 0;

}
int ConstructCovarianceMatrix(int NRate,int NFactor, double **C,int h)
{
	int i,j,k;

	for(i=0;i<NRate;i++)
	{
		for(j=0;j<=i;j++)
		{
			for(k=0 ; k<NFactor;k++)
			{
				C[i][j]+=A[i][k]*A[j][k];
			}
			C[j][i]=C[i][j];
		}

	}
	//AfficheMat(C,NRate,NRate);

	return 0;

}

double RecurseBushy(long h){
	long i,k;
	double tmp=0;
	double ProdLibor=1,SumLibor=0.0,SwapRate,B0=1,Bsi=1,SumBsi=0.0,res=0.0;
	int bushynoeud;
	SwapRate =0.0;
	bushynoeud=0;
	bushyprixnoeud = 0.0;
	NumeraireIndex = 0;
	if (h==NSteps){
		for(i=0;i<NRates;i++)
		{
			ProdLibor*=(1+EvolvedFra[h][i]*period);
			SumLibor += period/ProdLibor;
		}
		SwapRate=(1-(1/ProdLibor))/(SumLibor);

		for(i=0;i<NRates;i++)
		{
			Bsi*=1.0/(1+period*EvolvedFra[h][i]);
			SumBsi += Bsi*period;
		}
		res = (B0*SumBsi*MAX(type*(SwapRate-strike),0.0));
		return res;
	}

	for (i=0;i<NRates;i++){ // Calculate the drift for all rates and store them.
		mu_dT[h][i] = 0.;
		for (k=NumeraireIndex;k<=i;k++)
			mu_dT[h][i] += C[h][i][k] * EvolvedFra[h][k] * period / ( 1. + EvolvedFra[h][k] * (period) );
		for (k=i+1;k<NumeraireIndex;k++)
			mu_dT[h][i] -= C[h][i][k] * EvolvedFra[h][k] * period / ( 1. + EvolvedFra[h][k] * (period) );
	}
	if(h%bushytimeStep == 0)
	{
		for(i=0;i<NRates;i++)
		{
			ProdLibor*=(1+EvolvedFra[h][i]*period);
			SumLibor += period/ProdLibor;
		}
		SwapRate=(1-(1/ProdLibor))/(SumLibor);
		bushySwap[(long) (h/bushytimeStep)][bushyInd[(long) (h/bushytimeStep)]] = SwapRate;
	}

	for (k=0;k<NBranches;k++){ // Loop over all branches.

		for (i=0;i<NRates;i++){
			EvolvedFra[h+1][i] = EvolvedFra[h][i] * exp( mu_dT[h][i]*(Tau[h+1]-Tau[h]) + LogShiftOfBranch[h][k][i] );
		}

		if(h==0)
		{
			tmp+=RecurseBushy(h+1);
			h=0;
		}
		else if( (h%bushytimeStep == 0) && (h<=NSteps-bushytimeStep) && (bushytimeStep !=1) )
		{
			if(bushyInd[(long) (h/bushytimeStep)] % (long)( (pow(NBranches,bushytimeStep)-1)/(NBranches-1))== 0 || (bushynoeud==1) )
			{
				if(bushynoeud==0)
				{
					bushyInd[(long) (h/bushytimeStep)] +=1;
					bushyInd[(long)(h/bushytimeStep)+1]=0;
				}
				bushynoeud=1;
				tmp += RecurseBushy(h+1); // Sum up the results from all of the branches.
			}else
			{
				bushyInd[(long) (h/bushytimeStep)] +=1;
				return 0.0;
			}
		}else
		{
			if( (h>NSteps-bushytimeStep)|| bushytimeStep ==1)
			{
				tmp+=RecurseBushy(h+1);
				bushynoeud=(long)(bushytimeStep>1);
			}
			else
				RecurseBushy(h+1);
		}

	}

	if( (bushynoeud == 1 ) && ( h%bushytimeStep == 0 ))
	{
		int z=0;
		long tt= bushyInd[(long) (h/bushytimeStep)]/(long)( (pow(NBranches,bushytimeStep)-1)/(NBranches-1));
		if(bushyprixcalculer == 1)
		{
			bushyprix[(long) (h/bushytimeStep)][bushyInd[(long) (h/bushytimeStep)]-1] = bushyprixnoeud;
			bushyprixcalculer =0;
		}
		else
		{
			bushyprix[(long) (h/bushytimeStep)][bushyInd[(long) (h/bushytimeStep)]-1] = tmp;
		}

		if((tt == bushyBranches -1 )&&( bushytimeStep != 1 ) )
		{
			for(k=0;k<bushyBranches;k++)
			{
				x[k] = bushySwap[(long) (h/bushytimeStep)][(long)( k*( (pow(NBranches,bushytimeStep)-1)/(NBranches-1)) )];
				y[k] = bushyprix[(long) (h/bushytimeStep)][(long)(k*( (pow(NBranches,bushytimeStep)-1)/(NBranches-1)) )];
			}
			bushyprixnoeud=0.0;
			Sort(x,y,bushyBranches);
			for(k=0;k<(long)pow(NBranches,bushytimeStep);k++)
			{

				if(( k % (long) ((pow(NBranches,bushytimeStep)-1)/(NBranches-1)) )== 0)
					bushyprixnoeud+= y[(long)( k/(  (long)(pow(NBranches,bushytimeStep)-1)/(NBranches-1) ) )];
				else
					bushyprixnoeud+=cubspline(3,bushySwap[(long) (h/bushytimeStep)][k],x,y,bushyBranches);
			}
			bushyprixnoeud/=pow(NBranches,bushytimeStep);
			if(bushyprixnoeud<0)
				bushyprixnoeud=bushyprixnoeud;
			bushyprixcalculer = 1;
			return bushyprixnoeud;
		}else if(bushytimeStep == 1)
		{
			return tmp/NBranches;
		}
	}
	if( (h==0) && (k==NBranches))
	{
		if(bushytimeStep ==1)
			return tmp/NBranches;
		else
			return bushyprixnoeud/NBranches;
	}

	return tmp/NBranches;

}

double RecurseBushyCaplet(long h){
	long i,k;
	double tmp=0;
	double res=0.0;
	int bushynoeud;
	NumeraireIndex = 1;
	bushynoeud=0;
	bushyprixnoeud = 0.0;
	if (h==NSteps){
		res = (B0*period*MAX(type*(EvolvedFra[h][0]-strike),0.0));
		return res;
	}

	for (i=0;i<NRates;i++){ // Calculate the drift for all rates and store them.
		mu_dT[h][i] = 0.;
		for (k=NumeraireIndex;k<=i;k++)
			mu_dT[h][i] += C[h][i][k] * EvolvedFra[h][k] * period / ( 1. + EvolvedFra[h][k] * (period) );
		for (k=i+1;k<NumeraireIndex;k++)
			mu_dT[h][i] -= C[h][i][k] * EvolvedFra[h][k] * period / ( 1. + EvolvedFra[h][k] * (period) );
	}
	if(h%bushytimeStep == 0)
	{
		bushySwap[(long) (h/bushytimeStep)][bushyInd[(long) (h/bushytimeStep)]] = EvolvedFra[h][0];
	}

	for (k=0;k<NBranches;k++){ // Loop over all branches.

		for (i=0;i<NRates;i++){
			EvolvedFra[h+1][i] = EvolvedFra[h][i] * exp( mu_dT[h][i]*(Tau[h+1]-Tau[h]) + LogShiftOfBranch[h][k][i] );
		}

		if(h==0)
		{
			tmp+=RecurseBushyCaplet(h+1);
			h=0;
		}
		else if( (h%bushytimeStep == 0) && (h<=NSteps-bushytimeStep) && (bushytimeStep !=1) )
		{
			if(bushyInd[(long) (h/bushytimeStep)] % (long)( (pow(NBranches,bushytimeStep)-1)/(NBranches-1))== 0 || (bushynoeud==1) )
			{
				if(bushynoeud==0)
				{
					bushyInd[(long) (h/bushytimeStep)] +=1;
					bushyInd[(long)(h/bushytimeStep)+1]=0;
				}
				bushynoeud=1;
				tmp += RecurseBushyCaplet(h+1); // Sum up the results from all of the branches.
			}else
			{
				bushyInd[(long) (h/bushytimeStep)] +=1;
				return 0.0;
			}
		}else
		{
			if( (h>NSteps-bushytimeStep)|| bushytimeStep ==1)
			{
				tmp+=RecurseBushyCaplet(h+1);
				bushynoeud=(long)(bushytimeStep>1);
			}
			else
				RecurseBushyCaplet(h+1);
		}

	}

	if( (bushynoeud == 1 ) && ( h%bushytimeStep == 0 ))
	{
		int z=0;
		long tt= bushyInd[(long) (h/bushytimeStep)]/(long)( (pow(NBranches,bushytimeStep)-1)/(NBranches-1));
		if(bushyprixcalculer == 1)
		{
			bushyprix[(long) (h/bushytimeStep)][bushyInd[(long) (h/bushytimeStep)]-1] = bushyprixnoeud;
			bushyprixcalculer =0;
		}
		else
		{
			bushyprix[(long) (h/bushytimeStep)][bushyInd[(long) (h/bushytimeStep)]-1] = tmp;
		}
		if((tt == bushyBranches -1 )&&( bushytimeStep != 1 ) )
		{
			for(k=0;k<bushyBranches;k++)
			{
				x[k] = bushySwap[(long) (h/bushytimeStep)][(long)( k*( (pow(NBranches,bushytimeStep)-1)/(NBranches-1)) )];
				y[k] = bushyprix[(long) (h/bushytimeStep)][(long)(k*( (pow(NBranches,bushytimeStep)-1)/(NBranches-1)) )];
			}
			bushyprixnoeud=0.0;
			Sort(x,y,bushyBranches);
			for(k=0;k<(long)pow(NBranches,bushytimeStep);k++)
			{

				if(( k % (long) ((pow(NBranches,bushytimeStep)-1)/(NBranches-1)) )== 0)
					bushyprixnoeud+= y[(long)( k/(  (long)(pow(NBranches,bushytimeStep)-1)/(NBranches-1) ) )];
				else
					bushyprixnoeud+=cubspline(3,bushySwap[(long) (h/bushytimeStep)][k],x,y,bushyBranches);
			}
			bushyprixnoeud/=pow(NBranches,bushytimeStep);
			if(bushyprixnoeud<0)
				bushyprixnoeud=bushyprixnoeud;
			bushyprixcalculer = 1;
			return bushyprixnoeud;
		}else if(bushytimeStep == 1)
		{
			return tmp/NBranches;
		}
	}
	if( (h==0) && (k==NBranches))
	{
		if(bushytimeStep ==1)
			return tmp/NBranches;
		else
			return bushyprixnoeud/NBranches;
	}

	return tmp/NBranches;
}

double Recurse(long h){
	long i,k;
	double tmp=0;
	double ProdLibor=1,SumLibor=0.0,SwapRate,Bsi=1,SumBsi=0.0,res=0.0;
	NumeraireIndex= 0;
	if (h==NSteps){
		for(i=0;i<NRates;i++)
		{
			ProdLibor*=(1+EvolvedFra[h][i]*period);
			SumLibor += period/ProdLibor;
		}
		SwapRate=(1-(1/ProdLibor))/(SumLibor);

		for(i=0;i<NRates;i++)
		{
			Bsi*=1.0/(1+period*EvolvedFra[h][i]);
			SumBsi += Bsi;
		}
		res = (B0*SumBsi*period*MAX(type*(SwapRate-strike),0.0));
		return res;
		//return Intrinsic(h); // Termination of the recursion.
	}

	for (i=0;i<NRates;i++){ // Calculate the drift for all rates and store them.
		mu_dT[h][i] = 0.;
		for (k=NumeraireIndex;k<=i;k++)
			mu_dT[h][i] += C[h][i][k] * EvolvedFra[h][k] * period / ( 1. + EvolvedFra[h][k] * (period) );
	}
	for (k=0;k<NBranches;k++){ // Loop over all branches.
		for (i=0;i<NRates;i++){
			EvolvedFra[h+1][i] = EvolvedFra[h][i] * exp( mu_dT[h][i]*(Tau[h+1]-Tau[h]) + LogShiftOfBranch[h][k][i] );
		}
		tmp += Recurse(h+1); // Sum up the results from all of the branches.
	}

	return tmp/NBranches;
	//return CheckForEarlyExercise(h,tmp/NBranches);
}

double RecurseCaplet(long h){
	long i,k;
	double tmp=0;
	double res=0.0;
	NumeraireIndex = 1;
	if (h==NSteps){
		res = (B0*period*MAX(type*(EvolvedFra[h][0]-strike),0.0));
		return res;
	}

	for (i=0;i<NRates;i++){ // Calculate the drift for all rates and store them.
		mu_dT[h][i] = 0.;
		for (k=NumeraireIndex;k<=i;k++)
			mu_dT[h][i] += C[h][i][k] * EvolvedFra[h][k] * period / ( 1. + EvolvedFra[h][k] * (period) );
		for (k=i+1;k<NumeraireIndex;k++)
			mu_dT[h][i] -= C[h][i][k] * EvolvedFra[h][k] * period / ( 1. + EvolvedFra[h][k] * (period) );
	}
	for (k=0;k<NBranches;k++){ // Loop over all branches.
		for (i=0;i<NRates;i++){
			EvolvedFra[h+1][i] = EvolvedFra[h][i] * exp( mu_dT[h][i]*(Tau[h+1]-Tau[h]) + LogShiftOfBranch[h][k][i] );
		}

		tmp += RecurseCaplet(h+1); // Sum up the results from all of the branches.
	}
	return tmp/NBranches;
	//return CheckForEarlyExercise(h,tmp/NBranches);
}

double Interface()
{
	long NFactor,i,h,j,k;
	clock_t start, finish;
	double  duration;
	double price;
	NRates=4;
	NSteps=10;//nombre pas sur l'arbre et non le nombre de pas Ts/period comme dans MC
	period = 1;// c'est le Tau dans le libor exemple 0.5 pour 6 mois 1 pour 1 ans
	NBranches = 3; // c'est le nombre de branche sortant d'un noeud.
	NFactor=NBranches-1;// tjs comme ca à ne pas toucher

	bushytimeStep=2;// c'est le nombre de pas dans l'arbre unitaire de bushy non exlosif(tronqué)
	bushyBranches = NBranches;// c'est le nombre de branche  à partir
					//desquels on diffuse dans l'arbre bushy non exlosif(tronqué)
					// on peut aller jusqu'à pow(NBranches,bushytimeStep)
					//avec ces données on a 2 pas de temps dans Bushy unitaire cad on a à la fin de cet arbre
					// 9 noeud et on diffuse que à partir 3 noeuds : le plus bas  le plus haut et le milieu
					// sur les autres on interpole
	strike=0.075;
	Ts=4;
	type=-1;// 1 pour (S-K)+ et -1 pour (K-S)+
	typeOption = 1;// 1 pour Swaption et 0 pour (caplet ou florlet)
	flagBushy = 0;// si 0 methode de pricing sans decoupage de l'arbre;
			// si 1 methode de pricing avec un arbre non explosif decoupé
	flagRotate= 0; //0 : utilisation de la matrice de branchement S initiale
		// 1: rotation de la matrice de branchement S initiale
	B0=1;// ZeroCoupon
	// allocation de memoire et construction de la matrice de diffusion(de branchement de la gaussienne) avec rotation si flagRotate==1
	Allocate_All_memory();

	//mon libor initial;
	for(i=0;i<NRates;i++)
	{
		EvolvedFra[0][i] = 0.05;
	}
	// construction de la matrice de Covariance C apres la méthode de cholesky
	// C=A*A'
	// ATTENTION A LA DIFFERENCE ENTRE LE NSTEPS sur l'arbre et dans Monte-Carlo
	// si tu veux directement entrer la matrice decorrelée il faut directement remplir la matrice A
	// construction des termes martingales stockés dans LogShiftOfBranch en utilisant les termes gaussiens stockés dans Bbar
	for(h=0;h<NSteps;h++)
	{
		/*construire la matrice covariance pour une date dans **C1 et pour toutes les dates dans ***C ---*/
		//setCovarMatrix(C[h],h);
		/*decorreler pour chaque date h la matrice de covariance C[h], recuperer la matrice de cholesky A[][]*/
		//DecorrelateCovarianceMatrix(C[h],NRates,NFactor,C[h],h);//DecorrelateCovarianceMatrix(C[h],NRates,NFactor,C[h]);
		for(j=0;j<NRates;j++)
		{
			for( k=0;k<NFactor;k++)
			{
				if(k==0)
					A[j][k] = 0.2;
				else
				{
					A[j][k] = 0.1 + 0.07*((double)h/(NSteps-1));
				}
			}
		}
		ConstructCovarianceMatrix(NRates,NFactor,C[h],h);
		SetBbar(h);/*construire Bbar[][]*/
		SetLogShiftOfBranch(h);/*remplir le logshiftedOfBranch pour toute les dates*/
	}
	start = clock();
	if(typeOption == 0)
	{
		if(flagBushy ==0)
			price=Recurse(0);
		if(flagBushy ==1)
			price=RecurseBushy(0);
	}else
	{
		if(flagBushy ==0)
			price=RecurseCaplet(0);
		if(flagBushy ==1)
			price=RecurseBushyCaplet(0);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf( "%d \t %f \t %f\n",NSteps, duration,price );
	/*memory liberation*/
	Free_All_memory();
	return price;
}

double Price()
{
	double price=0.0;
	long NFactor,h=0,i;
	clock_t start, finish;
	double  duration;
	///bushytimeStep=2;
	//NSteps=16;
	//NRates=4;
	//NBranches=3;
	//NumeraireIndex=0;
	NFactor=NBranches-1;
	bushyBranches = NBranches;
	//strike=0.075;
	//Ts=4;
	//period=1;
	//type=-1;
	start = clock();

	/*memory allocation*/
	EvolvedFra_AllocatMem();
	C_AllocatMem();
	S_AllocatMem();
	A_AllocatMem();
	Tau_AllocatMem();
	mu_dT_AllocatMem();

	Bbar_AllocatMem();
	LogShiftOfBranch_AllocatMem();
	/*end memory allocation*/
	/*long *bushyInd;	  // de taille NSteps/bushytimeStep -1;
	double *bushySwap;// de taile NBranches^bushytimeStep
	double *bushyprix;// de taille NBranches^bushytimeStep*/
	bushyInd = (long *) malloc((long)(NSteps/bushytimeStep +1) * sizeof(long));
	for(h=0;h<(long)(NSteps/bushytimeStep +1) ;h++)
		bushyInd[h] = 0;
	bushySwap = (double **) malloc((long)(NSteps/bushytimeStep +1)*sizeof(double*));
	for(h=0;h<(long)(NSteps/bushytimeStep +1);h++)
	{
		bushySwap[h] = (double *) malloc((long)pow(NBranches,bushytimeStep)*sizeof(double));
		for(i=0;i<(long)pow(NBranches,bushytimeStep);i++)
			bushySwap[h][i] = 0.0;
	}

	bushyprix = (double **) malloc((long)(NSteps/bushytimeStep +1)*sizeof(double));
	for(h=0;h<(long)(NSteps/bushytimeStep +1) ;h++)
	{
		bushyprix[h] = (double *) malloc((long)pow(NBranches,bushytimeStep)*sizeof(double));
		for(i=0;i<(long)pow(NBranches,bushytimeStep);i++)
			bushyprix[h][i] = 0.0;
	}

	y= (double *) malloc((long)(bushyBranches)*sizeof(double));
	x= (double *) malloc((long)(bushyBranches)*sizeof(double));
	for(h=0;h<bushyBranches ;h++)
	{
		x[h] = 0.0;
		y[h] = 0.0;
	}

	SetS();
	//AfficheMat(S,NBranches,NBranches-1);
	if(flagRotate==1)
		RotateSimplex(NFactor);


	for(h=0;h<NSteps;h++)
	{
		/*construire la matrice covariance pour une date dans **C1 et pour toutes les dates dans ***C ---*/
		//setCovarMatrix(C[h],h);
		/*decorreler pour chaque date h la matrice de covariance C[h], recuperer la matrice de cholesky A[][]*/
		DecorrelateCovarianceMatrix(C[h],NRates,NFactor,C[h],h);//DecorrelateCovarianceMatrix(C[h],NRates,NFactor,C[h]);
		SetBbar(h);/*construire Bbar[][]*/
		SetLogShiftOfBranch(h);/*remplir le logshiftedOfBranch pour toute les dates*/
	}
	/*for(i=0;i<NRates;i++)
	{
		EvolvedFra[0][i] = 0.05;
	}*/
	EvolvedFra[0][0] = 0.06652;
	EvolvedFra[0][1] = 0.06251;
	EvolvedFra[0][2] = 0.06044;
	EvolvedFra[0][3] = 0.06044;
	/*appel a recurse	*/
	if(flagBushy==1)
		price=RecurseBushy(0);
	else
	{
		//price=RecurseBushyCaplet(0);
		//price=RecurseCaplet(0);
		price=Recurse(0);

	}
	/*memory liberation*/
	C_FreeMem();
	S_FreeMem();
	Tau_FreeMem();
	mu_dT_FreeMem();
	EvolvedFra_FreeMem();
	Bbar_FreeMem();
	LogShiftOfBranch_FreeMem();

	free(bushyInd);

	for(h=0;h<(long)(NSteps/bushytimeStep +1);h++)
	{
		free(bushySwap[h]);
		bushySwap[h]=NULL;
	}
	free(bushySwap);
	bushySwap=NULL;

	for(h=0;h<(long)(NSteps/bushytimeStep +1);h++)
	{
		free(bushyprix[h]);
		bushyprix[h]=NULL;
	}
	free(bushyprix);
	bushyprix=NULL;
	free(x);
	free(y);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf( "%d \t %f \t %f\n",NSteps, duration,price );
	fprintf(fic, "%d \t %f \t %f\n",NSteps, duration,price );
	return price;
}

void Cholesky(double **M, int dim)
{

	int i,j,k, erreur=0;
	double somme;
	i=0;
	while (!(erreur)&&(i<dim)){
		for (j=i;j<dim;j++){
			for (somme=M[i][j],k=i-1;k>=0;k--)
				somme-=M[i][k]*M[j][k];
			if (i==j) {
				if (somme<=0) erreur=1;
				else M[j][i]=sqrt(somme);
			}
			else M[j][i]=somme/M[i][i];
		}
		i++;
	}
	return;
}

int MyCholesky(double **M,double** B,int dim)
{
	long j,i,k;
	double somme=0.0;
	for( j= 0;j<dim;j++)
	{
		for(k=0;k<j-1;k++)
			somme +=B[j][k]*B[j][k];
		B[j][j]=sqrt(M[j][j]-somme);
		for(i=j;i<dim;i++)
		{
			somme=0;
			for(k=0;k<j-1;k++)
				somme+=B[i][k]*B[j][k];
			B[i][j]=(M[i][j]-somme)/B[j][j];
		}
	}

	return 0;
}

int AfficheMat(double **M,int diml,int dimc)
{
	long i,j;
	for(i=0;i<diml;i++)
	{
		for(j=0;j<dimc;j++)
		{
			printf("%f \t", M[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	return 0;

}
int main ()
{
	double price;
	price=Interface();
	return(1);
}
