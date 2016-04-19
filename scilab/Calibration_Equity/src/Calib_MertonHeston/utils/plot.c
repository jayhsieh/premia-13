#include <stdio.h>
#include <stdlib.h>
//#include "plot.h"

int plot(int sizeVector,  double *x, double *y)
{

  int i;
  FILE  *gp;
    

  gp = popen("gnuplot", "w");

  if(gp == NULL)
    {
      fprintf(stderr, "Oops, I can't find %s.", "gnuplot");
      exit(-1);
    }
  
  /*** initialization ***/
  fflush(gp);
  fprintf(gp, "set xlabel \"x\" \n");
  fprintf(gp, "set ylabel \"f(x)\" \n");
  //  fprintf(gp, "set pm3d\n");
  fprintf(gp, "plot '-' with point \n");
  
  for (i=0;i<sizeVector ;i++)
    {
	  fprintf(gp," %f %f \n",x[i],y[i]);
	  
	}

  fprintf(gp, "e \n");
  fflush(gp);
  printf("press_a_key_to_continue...\n");
  getchar();
  pclose(gp);

  return 0;
  
}



int plot3d(int Nx,int Ny,double *x,double *y,double **fxy)
{

  int i,j;
  FILE  *gp;
  FILE *fichier;
  

  gp = popen("gnuplot", "w");

  if(gp == NULL)
    {
      fprintf(stderr, "Oops, I can't find %s.", "gnuplot");
      exit(-1);
    }

  fichier = fopen ("surfprix1.dat", "wt");

  /*** initialization ***/
  fflush(gp);
  fprintf(gp, "set xlabel \"x\" \n");
  fprintf(gp, "set ylabel \"y\" \n");
  fprintf(gp, "set zlabel \"f(x,y)\" \n");
  fprintf(gp, "set hidden3d  \n");
  //  fprintf(gp, "set pm3d\n");
  fprintf(gp, "splot '-' with lines \n");
  
  for (i=0;i<Nx ;i++)
    {
	  for (j=0;j<Ny ;j++)
		{
		  fprintf(gp," %f %f %f \n",x[i],y[j],fxy[i][j]);
		  fprintf(fichier,"%f %f %f \n",x[i],y[j],fxy[i][j]);
		}
	  fprintf(gp," \n ");
	}
  

  fprintf(gp, "e \n");
  fflush(gp);
  printf("press_a_key_to_continue...\n");
  getchar();
  pclose(gp);
  fclose(fichier);
  
  
  return 0;
}

int plot3dbis(int Nx,int Ny,double *x,double *y,double *fxy)
{

  int i,j;
  FILE  *gp;
    

  gp = popen("gnuplot", "w");

  if(gp == NULL)
    {
      fprintf(stderr, "Oops, I can't find %s.", "gnuplot");
      exit(-1);
    }
  
  /*** initialization ***/
  fflush(gp);
  fprintf(gp, "set xlabel \"x\" \n");
  fprintf(gp, "set ylabel \"y\" \n");
  fprintf(gp, "set zlabel \"f(x,y)\" \n");
  //  fprintf(gp, "set pm3d\n");
  fprintf(gp, "splot '-' with lines \n");
  
  for (j=0;j<Ny ;j++)
    {
	  for (i=0;i<Nx ;i++)
		{
		  fprintf(gp," %f %f %f \n",x[i],y[j],fxy[j*Ny+i]);
		}
	  fprintf(gp," \n ");
	}
  

  fprintf(gp, "e \n");
  fflush(gp);
  printf("press_a_key_to_continue...\n");
  getchar();
  pclose(gp);
  
  return 0;
}

int myplot3d(int Nx,int Ny,double *x,double *y,double *fxy)
{

  int i,j;
  FILE  *gp;
  FILE *fichier;
  
  fichier = fopen("surfprix2.dat", "wt");
  
  gp = popen("gnuplot", "w");

  if(gp == NULL)
    {
      fprintf(stderr, "Oops, I can't find %s.", "gnuplot");
      exit(-1);
    }
  
  /*** initialization ***/
  fflush(gp);
  fprintf(gp, "set xlabel \"x\" \n");
  fprintf(gp, "set ylabel \"y\" \n");
  fprintf(gp, "set zlabel \"f(x,y)\" \n");
  //  fprintf(gp, "set pm3d\n");
  fprintf(gp, "splot '-' with lines \n");
  
  for (j=0;j<Ny ;j++)
    {
	  for (i=0;i<Nx ;i++)
		{
		  fprintf(gp," %f %f %f \n",x[i],y[j],fxy[j*Ny+i+1]);
		  fprintf(fichier," %f %f %f \n",x[i],y[j],fxy[j*Ny+i+1]);
		}
	  fprintf(gp," \n ");
	}
  

  fprintf(gp, "e \n");
  fflush(gp);
  printf("press_a_key_to_continue...\n");
  getchar();
  pclose(gp);
  fclose(fichier);
 
  return 0;
}

