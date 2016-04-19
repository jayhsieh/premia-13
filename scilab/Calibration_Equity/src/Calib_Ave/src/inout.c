#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "inout.h"


/*
  MATHFI Project, Inria Rocquencourt.

*/

void loadParamPricer(char *name_in_pricer,double *pt_K, double *pt_T, char **pt_option_type,char ** pt_name_in_grid_price_in,char
**pt_name_in_grid_price_out)
{
  char *ligne;
  char *ok;

  int cpt;
  double *param_double[2];
  char **param_char[3];

  int cpt_double;
  int cpt_char;
  int nbr;
  int i=0;

  FILE *fic_in_pricer;
  fic_in_pricer = fopen(name_in_pricer,"r");

  param_double[0] = pt_K;
  param_double[1] = pt_T;

  param_char[0] = pt_option_type;
  param_char[1] = pt_name_in_grid_price_in;
  param_char[2] = pt_name_in_grid_price_out;
  cpt_double = 0;
  cpt_char = 0;

  nbr = 100;
  //3 comme deux parametre char**
  for (i=0;i<3;i++)
    *(param_char[i]) = (char *) malloc(nbr*sizeof(char));

  ligne = (char *) malloc(nbr*sizeof(char));
	// 5 parametres
  for (cpt=1;cpt<=5;cpt++){

    ok = fgets(ligne,nbr,fic_in_pricer);
    if (ok == NULL){
      printf("Pb de lecture dans le fichier d'entree\n");
      exit(-1);
    }

    while (ligne[0] == '#'){
      ok = fgets(ligne,nbr,fic_in_pricer);
      if (ok == NULL){
	printf("Pb de lecture dans le fichier d'entree\n");
	exit(-1);
      }
    }

    /* dans ligne : la valeur numero cpt+1 sous forme de chaine de carac */
    if (cpt==3 || cpt==4 || cpt==5){
      /* on doit convertir la chaine en string */
      stringToString(*(param_char[cpt_char]),ligne);
      cpt_char++;
    }
    else{
      /* on doit convertir la chaine en double */
      *(param_double[cpt_double]) = stringToDouble(ligne);
      cpt_double++;
    }

  }

  fclose(fic_in_pricer);

  free(ligne);


}

void loadParamCalib(char *name_in_calib, double *pt_S_0, double *pt_r, double *pt_q, int *pt_N ,\
 double *pt_sigma_0, double *pt_sigma_min, double *pt_sigma_max, double *pt_sigma_bar, \
 double *pt_gradtol,double *pt_steptol,int *pt_verbosity,int *pt_saveSuccessiveXinFile,\
 int *pt_maxCounter,double *pt_lambda, double *pt_alpha,char **pt_name_in_data,char **pt_name_out_Vol_Loc)
{

  char *ligne;
  char *ok;
  int nbr;
  int cpt;
  double *param_double[11];
  int *param_int[4];
  char **param_char[2];

  int cpt_double;
  int cpt_int;
  int cpt_char;

  int i=0;

  FILE *fic_in_calib;
  fic_in_calib = fopen(name_in_calib,"r");

  param_double[0] = pt_S_0;
  param_double[1] = pt_r;
  param_double[2] = pt_q;
  param_double[3] = pt_sigma_0;
  param_double[4] = pt_sigma_min;
  param_double[5] = pt_sigma_max;
  param_double[6] = pt_sigma_bar;
  param_double[7] = pt_gradtol ;
  param_double[8] = pt_steptol ;
  param_double[9] = pt_lambda ;
  param_double[10] = pt_alpha ;

  param_int[0] = pt_N;
  param_int[1] = pt_verbosity;
  param_int[2] = pt_saveSuccessiveXinFile;
  param_int[3] = pt_maxCounter;

  param_char[0] = pt_name_in_data;
  param_char[1] = pt_name_out_Vol_Loc;

  cpt_double = 0;
  cpt_int = 0;
  cpt_char = 0;

  nbr = 100;
  //2 comme deux parametre char**
  for (i=0;i<2;i++)
    *(param_char[i]) = (char *) malloc(nbr*sizeof(char));

  ligne = (char *) malloc(nbr*sizeof(char));
	// 17 parametres
  for (cpt=1;cpt<=17;cpt++){

    ok = fgets(ligne,nbr,fic_in_calib);
    if (ok == NULL){
      printf("Pb de lecture dans le fichier d'entree\n");
      exit(-1);
    }

    while (ligne[0] == '#'){
      ok = fgets(ligne,nbr,fic_in_calib);
      if (ok == NULL){
	printf("Pb de lecture dans le fichier d'entree\n");
	exit(-1);
      }
    }

    /* dans ligne : la valeur numero cpt+1 sous forme de chaine de carac */
    if (cpt==4||cpt==11||cpt==12||cpt==13){
      /* on doit convertir la chaine en entier */
      *(param_int[cpt_int]) = stringToInt(ligne);
      cpt_int++;
    }
    else if (cpt==16 || cpt == 17){
      stringToString(*(param_char[cpt_char]),ligne);
      cpt_char++;
    }

    else{
      /* on doit convertir la chaine en double */
      *(param_double[cpt_double]) = stringToDouble(ligne);
      cpt_double++;
    }

  }

  fclose(fic_in_calib);

  free(ligne);

}


void stringToString(char *strcopy, char *string){

  sscanf(string,"%s",strcopy);

}

int stringToInt(char *string){

  /*string se termine par \n\0 */
  int res=0;
  sscanf(string,"%i",&res);

  return res;

}

double stringToDouble(char *string){

  double res=0;
  sscanf(string,"%lf",&res);

    return res;

}

