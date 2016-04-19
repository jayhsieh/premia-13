#include "MMsimulBGM.h"

int C_AllocatMem(){
	long h,i,j;
	C= (double ***) malloc(NSteps*sizeof(double**));
	for (h=0;h<NSteps;h++){
		C[h] =(double **) malloc(NRates*sizeof(double*));
		for (i=0;i<NRates;i++){
			C[h][i] =(double*) malloc(NRates*sizeof(double));
			for(j=0;j<NRates;j++)
				C[h][i][j]=0.0;
		}
	}
	return 0;
	
}
  //mu_dT[i] double * de taille NRates; 
  //EvolvedFra[h][k] de taille [NSteps*NRates]
int EvolvedFra_AllocatMem(){
	long h;
	EvolvedFra=(double **) malloc((NSteps+1)*sizeof(double));
	for (h=0;h<NSteps+1;h++){
		EvolvedFra[h] =(double*) malloc(NRates*sizeof(double));
	}
	return 0;
}
int A_AllocatMem(){
	long h;
	A=(double **) malloc(NRates*sizeof(double));
	for (h=0;h<NRates;h++){
		A[h] =(double*) malloc((NBranches-1)*sizeof(double));
	}
	return 0;
}

int Bbar_AllocatMem(){
	long h,i,k;
	Bbar =(double ***) malloc(NSteps*sizeof(double));
	for(h=0;h<NSteps;h++)
	{
		Bbar[h]=(double **) malloc(NRates*sizeof(double));
		for (i=0;i<NRates;i++){
			Bbar[h][i] =(double*) malloc(NBranches*sizeof(double));
			for(k=0;k<NBranches;k++)
				Bbar[h][i][k]=0.0;

		}
	}
	return 0;
}
int S_AllocatMem(){
	long h;
	S=(double **) malloc(NBranches*sizeof(double));
	for (h=0;h<NBranches;h++){
		S[h] =(double*) malloc((NBranches-1)*sizeof(double));
	}
	return 0;
}

/*------------------Tau[k] de taille[NRates]; c une approximation grossaire ---------------------*/
/* ----------normalement pour les rabres il faut rafiner le pas de temps-------------------------*/
int Tau_AllocatMem(){
	int i;
	Tau =(double *) malloc((NSteps+1)*sizeof(double));
	Tau[0]=0.0;
	for(i=1;i<NSteps+1;i++)
	{
		Tau[i] =Tau[i-1] + Ts/((double)NSteps);
	}
	return 0;
}
int mu_dT_AllocatMem(){
	long h;
	mu_dT =(double **) malloc(NSteps*sizeof(double));
	for (h=0;h<NSteps;h++){
		mu_dT[h] =(double *) malloc(NRates*sizeof(double));
	}
	return 0;
}
//LogShiftOfBranch[h][k][i] de taille [NSteps][NRatess][NRates];
int LogShiftOfBranch_AllocatMem(){
	long h,k;
	LogShiftOfBranch= (double ***) malloc(NSteps*sizeof(double**));
	for (h=0;h<NSteps;h++){
		LogShiftOfBranch[h] =(double **) malloc(NBranches*sizeof(double*));
		for (k=0;k<NBranches;k++){
			LogShiftOfBranch[h][k] =(double*) malloc(NRates*sizeof(double));
		}
	}
	return 0;
}

int C_FreeMem()
{
	long h,i;
	for (h=0;h<NSteps;h++){
		for (i=0;i<NRates;i++){
			free(C[h][i]);
			C[h][i]=NULL;
		}
		free(C[h]);
		C[h]=NULL;
	}
	free(C);
	C=NULL;
	
	return 0;
}
int A_FreeMem(){
	long h;
	
	for (h=0;h<NRates;h++){
		free(A[h]);
		A[h]=NULL;
	}
	free(A);
	A=NULL;
	return 0;
}

int	S_FreeMem()
{
	long h;
	for (h=0;h<NBranches;h++){
		free(S[h]);
		S[h]=NULL;
	}
	free(S);
	S=NULL;
	return 0;
}

int Tau_FreeMem()
{
	free(Tau);
	Tau=NULL;
	return 0;
}
int mu_dT_FreeMem()
{
	long h;
	for (h=0;h<NSteps;h++){
		free(mu_dT[h]);
		mu_dT[h]=NULL;
	}
	free(mu_dT);
	mu_dT=NULL;
	return 0;
}
int	EvolvedFra_FreeMem()
{
	long h;
	for (h=0;h<NSteps+1;h++){
		free(EvolvedFra[h]);
		EvolvedFra[h]=NULL;
	}
	free(EvolvedFra);
	EvolvedFra=NULL;
	return 0;
}
int	Bbar_FreeMem()
{
	long h,i;
	for(h=0;h<NSteps;h++)
	{
		for (i=0;i<NRates;i++){
			free(Bbar[h][i]);
			Bbar[h][i]=NULL;
		}
		free(Bbar[h]);
		Bbar[h]=NULL;
	}
	free(Bbar);
	Bbar=NULL;
	return 0;
}
int	LogShiftOfBranch_FreeMem()
{
	long h,k;
	for (h=0;h<NSteps;h++){
		for (k=0;k<NBranches;k++){
			free(LogShiftOfBranch[h][k]);
			LogShiftOfBranch[h][k]=NULL;
		}
		free(LogShiftOfBranch[h]);
		LogShiftOfBranch[h]=NULL;
	}
	free(LogShiftOfBranch);
	LogShiftOfBranch=NULL;
	return 0;

}
int Allocate_All_memory()
{
	long h,i;
	EvolvedFra_AllocatMem();
	C_AllocatMem();
	S_AllocatMem();
	A_AllocatMem();
	Tau_AllocatMem();
	mu_dT_AllocatMem();

	Bbar_AllocatMem();
	LogShiftOfBranch_AllocatMem();
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
	if(flagRotate==1)
		RotateSimplex(NBranches-1);
	return 0;
}
int Free_All_memory()
{
	long h,i;
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
	return 0;
}


