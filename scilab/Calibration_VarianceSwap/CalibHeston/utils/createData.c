 
#include <stdio.h>
#include <stdlib.h>


int main(int argc , char ** argv)
{

char nom;
int i, nbstrikes, nbmat, type, curr, firstT, lastT, typen, nbetapes, nbdata;
double T, K, sigma, r, d, wi, St0;
double sigma0, V0, kappa, theta, sigmav, rho, fichier, matur;
double *maturs, *strikes, **ivol;
FILE *fich;

	nbmat=14;
	nbstrikes=21;
	printf("Number of maturities: %d\n Number of strikes: %d\n", nbmat, nbstrikes);

	maturs = (double*)malloc(sizeof(double)*nbmat);
	strikes = (double*)malloc(sizeof(double)*nbstrikes);

	ivol = (double **)malloc(nbmat*sizeof(double *));
	for (i=0;i<nbmat;i++)	ivol[i] = (double *) malloc(nbstrikes*sizeof(double));

	printf("Reading file: PrixOption.txt...");
//------------------------------PRIXOPTION.TXT------------------------------
	fich = fopen ("PrixOption.txt", "r");
	for(curr=0; curr<nbmat; curr++)
	{
		fscanf(fich, "%lf", &maturs[curr]);
	}
	for(i=0;i<nbstrikes;i++)
	{
		fscanf(fich, "%lf", &strikes[i]);
		for(curr=0; curr<nbmat; curr++)
		{
			fscanf(fich, "%lf", &ivol[curr][i]);
		}
	}
	fclose(fich);
//--------------------------------------------------------------------------
	printf("Done\n");
	firstT=6;
	lastT=11;
	
	printf("Write to DataMarket.txt...");
//------------------------------DATAMARKET.TXT------------------------------
	fich = fopen ("DataMarket.txt", "w");
	for(curr=firstT; curr<lastT+1; curr++)
	{
		T=maturs[curr];
		for(i=0;i<nbstrikes;i++)
		{
	  	//
			K=strikes[i];
			type=K<1? -1: 1;
			sigma=ivol[curr][i];
			r=0.0;
			d=0.0;
			wi=1.0;
	  		fprintf(fich,"%d %lf %lf %lf %lf %lf %lf\n", type, T, K, sigma, r, d, wi);
		}
	}
	fclose(fich);
//--------------------------------------------------------------------------
	printf("Done\n");

	typen=33;
	nbetapes=1;
	nbdata=nbstrikes*(lastT-firstT+1);
	sigma0  = 0.3;
	V0      = sigma0*sigma0;
	kappa   = 2.; 
    	theta   = 0.04; 
    	sigmav  = 0.3;
    	rho     = -0.5;
	
	printf("Write to in.dat...");
//-----------------------------------------------IN.DAT----------------------	
	fich = fopen ("in.dat","w");
  //
	fprintf(fich,"%d \n",typen);
  	fprintf(fich,"%d \n",nbetapes);
  	fprintf(fich,"%d \n",nbdata);
 	fprintf(fich,"%lf %lf %lf %lf %lf\n", V0, kappa, theta, sigmav, rho);
    
	fclose(fich);
//------------------------------------------------------------------------------
	printf("Done\n");
	
	for (i=0;i<nbmat;i++)	free(ivol[i]);
	free(ivol);
	free(maturs);
	free(strikes);

  //
  return 0;
  //
}
