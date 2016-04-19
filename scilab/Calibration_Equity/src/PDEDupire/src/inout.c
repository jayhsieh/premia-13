#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "inout.h"
#include "spline.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle and Jean-Marc Cognet, November 2002.
*/

void loadSigmaddl(char *name_in_sigma, int *pt_n, int *pt_m, double **pt_y_coarseGrid, double **pt_T_coarseGrid, double **pt_sigma){

  int nbParamSigma;
  int i;
  FILE *fic_in_sigma;

  fic_in_sigma = fopen(name_in_sigma,"r"); 

  fscanf(fic_in_sigma,"%d",pt_n);
  fscanf(fic_in_sigma,"%d",pt_m);
  nbParamSigma = (*pt_n+3)*(*pt_m+3);

  *pt_y_coarseGrid = (double *)malloc((*pt_n+1)*sizeof(double));
  *pt_T_coarseGrid = (double *)malloc((*pt_m+1)*sizeof(double));
  *pt_sigma = (double *)malloc(nbParamSigma*sizeof(double));

  for (i=0;i<*pt_n+1;i++)
    fscanf(fic_in_sigma,"%lf",&((*pt_y_coarseGrid)[i]));

  for (i=0;i<*pt_m+1;i++)
    fscanf(fic_in_sigma,"%lf",&((*pt_T_coarseGrid)[i]));
  
  for (i=0;i<nbParamSigma;i++) 
    fscanf(fic_in_sigma,"%lf",&((*pt_sigma)[i]));

  fclose(fic_in_sigma); 

}


void loadSigmaGrid(char *name_in_sigma, int *pt_n, int *pt_m, double **pt_k_coarseGrid, double **pt_T_coarseGrid){

  int i;
  FILE *fic_in_sigma;

  fic_in_sigma = fopen(name_in_sigma,"r"); 

  fscanf(fic_in_sigma,"%d",pt_n);
  fscanf(fic_in_sigma,"%d",pt_m);
 
  *pt_k_coarseGrid = (double *)malloc((*pt_n+1)*sizeof(double));
  *pt_T_coarseGrid = (double *)malloc((*pt_m+1)*sizeof(double));
 
  for (i=0;i<*pt_n+1;i++)
    fscanf(fic_in_sigma,"%lf",&((*pt_k_coarseGrid)[i]));

  for (i=0;i<*pt_m+1;i++)
    fscanf(fic_in_sigma,"%lf",&((*pt_T_coarseGrid)[i]));
  
  fclose(fic_in_sigma); 

}

void sigmaddlToSigmaFineGrid(double **sigmaFineGrid, int n, int m, double *y_coarseGrid, double *T_coarseGrid, double *sigma_param, double *y_fineGrid, double *T_fineGrid, int N, int M){

  int i,k,l;
  double **sigmaCoarseGrid;
  struct derivData *interpolData;

  sigmaCoarseGrid = (double **) malloc((n+1)*sizeof(double *));
  for (k=0;k<n+1;k++){
    sigmaCoarseGrid[k] = (double *) malloc((m+1)*sizeof(double));
    for (l=0;l<m+1;l++)
      sigmaCoarseGrid[k][l] = sigma_param[l*(n+1)+k];
  }

  interpolData = (struct derivData *) malloc(sizeof(struct derivData));
  interpolData->deriv_y_0 = (double *) malloc((m+1)*sizeof(double));
  interpolData->deriv_y_n = (double *) malloc((m+1)*sizeof(double));
  interpolData->deriv_T_0 = (double *) malloc((n+1)*sizeof(double));
  interpolData->deriv_T_m = (double *) malloc((n+1)*sizeof(double));

  for (k=(n+1)*(m+1);k<(n+1)*(m+1)+m+1;k++)
    interpolData->deriv_y_0[k-(n+1)*(m+1)] = sigma_param[k];

  for (k=(n+2)*(m+1);k<(n+2)*(m+1)+m+1;k++)
    interpolData->deriv_y_n[k-(n+2)*(m+1)] = sigma_param[k];
  
  for (k=(n+3)*(m+1);k<(n+3)*(m+1)+n+1;k++)
    interpolData->deriv_T_0[k-(n+3)*(m+1)] = sigma_param[k];

  for (k=(n+3)*(m+1)+n+1;k<(n+3)*(m+1)+2*(n+1);k++)
    interpolData->deriv_T_m[k-((n+3)*(m+1)+n+1)] = sigma_param[k];

  interpolData->deriv_yT_00 = sigma_param[(n+3)*(m+1)+2*(n+1)];
  interpolData->deriv_yT_n0 = sigma_param[(n+3)*(m+1)+2*(n+1)+1];
  interpolData->deriv_yT_0m = sigma_param[(n+3)*(m+1)+2*(n+1)+2];
  interpolData->deriv_yT_nm = sigma_param[(n+3)*(m+1)+2*(n+1)+3];

  interpolData->n = n;
  interpolData->m = m;

  /* interpolation of sigmaCoarseGrid to sigmaFineGrid */
  
  interpole(sigmaFineGrid,sigmaCoarseGrid,n,m,N,M,interpolData,y_coarseGrid,T_coarseGrid,y_fineGrid,T_fineGrid);

}

void loadDataPrices(char *name_in_data, struct marketData **data){

  FILE *fic_in_data;
  double K,T,V;
  struct marketData *pt_data;
  
  fic_in_data = fopen(name_in_data,"r");

  if (fscanf(fic_in_data,"%lf %lf %lf",&K,&T,&V) == 3){
    (*data) = (struct marketData *) malloc(sizeof(struct marketData)); 
    (*data)->strike = K;
    (*data)->maturity = T;
    (*data)->price = V;
    (*data)->next = NULL;
    pt_data = (*data);

    while (fscanf(fic_in_data,"%lf %lf %lf",&K,&T,&V) == 3){
      pt_data->next = (struct marketData *) malloc(sizeof(struct marketData));
      pt_data = pt_data->next;
      pt_data->strike = K;
      pt_data->maturity = T;
      pt_data->price = V;
      pt_data->next = NULL;
    }  
  
  }
  else
    (*data) = NULL;
      
   fclose(fic_in_data);

}

//void loadSigmaddl(char *name_in_sigma_ddl, int *pt_n, int *pt_m, double **y_coarseGrid, double **T_coarseGrid, double **sigma_ddl){
/* void loadSigmaddl(char *name_in_sigma_ddl, int *pt_n, int *pt_m, double **pt_y_coarseGrid){ */

/*   char *ligne; */
/*   char *ok; */
/*   int i, nbr; */

/*   FILE *fic_in_sigma_ddl; */

/*   nbr = 100; */
/*   ligne = (char *) malloc(nbr*sizeof(char)); */

/*   fic_in_sigma_ddl = fopen(name_in_sigma_ddl,"r"); */

/*   ok = fgets(ligne,nbr,fic_in_sigma_ddl); */
/*   *pt_n = stringToInt(ligne); */
/*   ok = fgets(ligne,nbr,fic_in_sigma_ddl); */
/*   *pt_m = stringToInt(ligne); */

/*   (*pt_y_coarseGrid) = (double *) malloc((*pt_n+1)*sizeof(double)); */

/*   for (i=0;i<*pt_n+1;i++){ */
/*     ok = fgets(ligne,nbr,fic_in_sigma_ddl); */
/*     *(pt_y_coarseGrid[i]) = stringToDouble(ligne); */
/*   } */

/*   fclose(fic_in_sigma_ddl); */

/*   free(ligne); */

/* } */

void loadParamSimul(char *name_in_simul, double *pt_S_0, double *pt_r, double *pt_q, int *pt_optionType, int *pt_optionSimul, double *pt_t_0, double *pt_T_max, double *pt_y_min, double *pt_y_max, int *pt_N, int *pt_M, int *pt_gridType, double *pt_theta, char **pt_name_in_sigma_ddl, double *pt_sigmaCte, char **pt_name_in_visu, char **pt_name_out_visu, char **pt_name_in_data, char **pt_name_out_data){

  char *ligne;
  char *ok;
  int nbr;
  int cpt;
  double *param_double[9];
  int *param_int[5];  
  char **param_char[5];

  int cpt_double;
  int cpt_int;
  int cpt_char;

  int i=0;

  FILE *fic_in_simul;
  fic_in_simul = fopen(name_in_simul,"r");

  param_double[0] = pt_S_0;
  param_double[1] = pt_r;
  param_double[2] = pt_q;
  param_double[3] = pt_t_0;
  param_double[4] = pt_T_max;
  param_double[5] = pt_y_min;
  param_double[6] = pt_y_max;
  param_double[7] = pt_theta;
  param_double[8] = pt_sigmaCte;

  param_int[0] = pt_optionType;
  param_int[1] = pt_optionSimul;
  param_int[2] = pt_N;
  param_int[3] = pt_M;
  param_int[4] = pt_gridType;
  
  param_char[0] = pt_name_in_sigma_ddl;
  param_char[1] = pt_name_in_visu;
  param_char[2] = pt_name_out_visu;
  param_char[3] = pt_name_in_data;
  param_char[4] = pt_name_out_data;

  cpt_double = 0;
  cpt_int = 0;
  cpt_char = 0;

  nbr = 100;
  
  for (i=0;i<5;i++)
    *(param_char[i]) = (char *) malloc(nbr*sizeof(char));
  
  ligne = (char *) malloc(nbr*sizeof(char));

  for (cpt=1;cpt<=19;cpt++){

    ok = fgets(ligne,nbr,fic_in_simul);
    if (ok == NULL){
      printf("Pb de lecture dans le fichier d'entree\n");
      exit;
    }
    
    while (ligne[0] == '#'){
      ok = fgets(ligne,nbr,fic_in_simul);
      if (ok == NULL){
	printf("Pb de lecture dans le fichier d'entree\n");
	exit;
      }
    }
    
    /* dans ligne : la valeur numero cpt+1 sous forme de chaine de carac */
    if (cpt==4 || cpt==5 || cpt==10 || cpt==11 || cpt==12){
      /* on doit convertir la chaine en entier */
      *(param_int[cpt_int]) = stringToInt(ligne);
      cpt_int++;
    }
    else if (cpt==14 || cpt==16 || cpt==17 || cpt==18 || cpt==19){ 
      stringToString(*(param_char[cpt_char]),ligne);
      cpt_char++;
    }
    else{
      /* on doit convertir la chaine en double */      
      *(param_double[cpt_double]) = stringToDouble(ligne);      
      cpt_double++;
    }

  }

  fclose(fic_in_simul);

  free(ligne);

}

void loadParamCalib(char *name_in_calib, double *pt_S_0, double *pt_r, double *pt_q, int *pt_optionType, double *pt_t_0, double *pt_T_max, double *pt_y_min, double *pt_y_max, int *pt_N, int *pt_M, int *pt_gridType, double *pt_theta, int *pt_choice_optim, char **pt_name_in_optim, char **pt_name_in_data, char **pt_name_in_sigmainit_ddl, char **pt_name_out_sigmaest_ddl, char **pt_name_in_sigma_visu, char **pt_name_out_sigmainit_visu, char **pt_name_out_sigmaest_visu){

  char *ligne;
  char *ok;
  int nbr;
  int cpt;
  double *param_double[8];
  int *param_int[5];
  char **param_char[7];

  int cpt_double;
  int cpt_int;
  int cpt_char;

  int i=0;

  FILE *fic_in_calib;
  fic_in_calib = fopen(name_in_calib,"r");

  param_double[0] = pt_S_0;
  param_double[1] = pt_r;
  param_double[2] = pt_q;
  param_double[3] = pt_t_0;
  param_double[4] = pt_T_max;
  param_double[5] = pt_y_min;
  param_double[6] = pt_y_max;
  param_double[7] = pt_theta;

  param_int[0] = pt_optionType;
  param_int[1] = pt_N;
  param_int[2] = pt_M;
  param_int[3] = pt_gridType;
  param_int[4] = pt_choice_optim;
  
  param_char[0] = pt_name_in_optim;
  param_char[1] = pt_name_in_data;
  param_char[2] = pt_name_in_sigmainit_ddl;
  param_char[3] = pt_name_out_sigmaest_ddl;
  param_char[4] = pt_name_in_sigma_visu;
  param_char[5] = pt_name_out_sigmainit_visu;
  param_char[6] = pt_name_out_sigmaest_visu;

  cpt_double = 0;
  cpt_int = 0;
  cpt_char = 0;

  nbr = 100;
  
  for (i=0;i<7;i++)
    *(param_char[i]) = (char *) malloc(nbr*sizeof(char));
  
  ligne = (char *) malloc(nbr*sizeof(char));

  for (cpt=1;cpt<=20;cpt++){

    ok = fgets(ligne,nbr,fic_in_calib);
    if (ok == NULL){
      printf("Pb de lecture dans le fichier d'entree\n");
      exit;
    }
    
    while (ligne[0] == '#'){
      ok = fgets(ligne,nbr,fic_in_calib);
      if (ok == NULL){
	printf("Pb de lecture dans le fichier d'entree\n");
	exit;
      }
    }
    
    /* dans ligne : la valeur numero cpt+1 sous forme de chaine de carac */
    if (cpt==4 || cpt==9 || cpt==10 || cpt==11 || cpt==13){
      /* on doit convertir la chaine en entier */
      *(param_int[cpt_int]) = stringToInt(ligne);
      cpt_int++;
    }
    else if (cpt>=14){ 
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

void loadParamRafsig(char *name_in_rafsig, char **pt_name_in_sigma1_ddl, char **pt_nbsplit_y_char, char **pt_nbsplit_T_char, char **pt_name_out_sigma2_ddl){

  char *ligne;
  char *ok;
  int nbr;
  int cpt;
  char **param_char[4];

  int i=0;

  FILE *fic_in_rafsig;
  fic_in_rafsig = fopen(name_in_rafsig,"r");

  param_char[0] = pt_name_in_sigma1_ddl;
  param_char[1] = pt_nbsplit_y_char;
  param_char[2] = pt_nbsplit_T_char;
  param_char[3] = pt_name_out_sigma2_ddl;

  nbr = 100;
  
  for (i=0;i<4;i++)
    *(param_char[i]) = (char *) malloc(nbr*sizeof(char));
  
  ligne = (char *) malloc(nbr*sizeof(char));

  for (cpt=0;cpt<4;cpt++){

    ok = fgets(ligne,nbr,fic_in_rafsig);
    if (ok == NULL){
      printf("Pb de lecture dans le fichier d'entree\n");
      exit;
    }
    
    while (ligne[0] == '#'){
      ok = fgets(ligne,nbr,fic_in_rafsig);
      if (ok == NULL){
	printf("Pb de lecture dans le fichier d'entree\n");
	exit;
      }
    }
    
    stringToString( *(param_char[cpt]),ligne); 
    if (cpt==1 || cpt==2)
      strcat(*(param_char[cpt]),"\n");

  }

  fclose(fic_in_rafsig);

  free(ligne);

}

void loadParamVisusig(char *name_in_visusig, char **pt_name_in_sigma_ddl, char **pt_name_in_sigma_visu, char **pt_name_out_sigma_visu){

  char *ligne;
  char *ok;
  int nbr;
  int cpt;
  char **param_char[3];

  int i=0;

  FILE *fic_in_visusig;
  fic_in_visusig = fopen(name_in_visusig,"r");

  param_char[0] = pt_name_in_sigma_ddl;
  param_char[1] = pt_name_in_sigma_visu;
  param_char[2] = pt_name_out_sigma_visu;

  nbr = 100;
  
  for (i=0;i<3;i++)
    *(param_char[i]) = (char *) malloc(nbr*sizeof(char));
  
  ligne = (char *) malloc(nbr*sizeof(char));

  for (cpt=0;cpt<3;cpt++){

    ok = fgets(ligne,nbr,fic_in_visusig);
    if (ok == NULL){
      printf("Pb de lecture dans le fichier d'entree\n");
      exit;
    }
    
    while (ligne[0] == '#'){
      ok = fgets(ligne,nbr,fic_in_visusig);
      if (ok == NULL){
	printf("Pb de lecture dans le fichier d'entree\n");
	exit;
      }
    }
    
    stringToString( *(param_char[cpt]),ligne); 

  }

  fclose(fic_in_visusig);

  free(ligne);

}

void loadParamImpsig(char *name_in_impsig, double *pt_S_0, double *pt_r, double *pt_q, int *pt_optionType, double *pt_t_0, char **pt_name_in_data, char **pt_name_out_sigma){

  char *ligne;
  char *ok;
  int nbr;
  int cpt;
  double *param_double[4];
  int *param_int;  
  char **param_char[2];

  int cpt_double;
  int cpt_char;

  int i=0;

  FILE *fic_in_impsig;
  fic_in_impsig = fopen(name_in_impsig,"r");

  param_double[0] = pt_S_0;
  param_double[1] = pt_r;
  param_double[2] = pt_q;
  param_double[3] = pt_t_0;

  param_int = pt_optionType;
  
  param_char[0] = pt_name_in_data;
  param_char[1] = pt_name_out_sigma;


  cpt_double = 0;
  cpt_char = 0;

  nbr = 100;
  
  for (i=0;i<2;i++)
    *(param_char[i]) = (char *) malloc(nbr*sizeof(char));
  
  ligne = (char *) malloc(nbr*sizeof(char));

  for (cpt=1;cpt<=7;cpt++){

    ok = fgets(ligne,nbr,fic_in_impsig);
    if (ok == NULL){
      printf("Pb de lecture dans le fichier d'entree\n");
      exit;
    }
    
    while (ligne[0] == '#'){
      ok = fgets(ligne,nbr,fic_in_impsig);
      if (ok == NULL){
	printf("Pb de lecture dans le fichier d'entree\n");
	exit;
      }
    }
    
    /* dans ligne : la valeur numero cpt+1 sous forme de chaine de carac */
    if (cpt==4)
      /* on doit convertir la chaine en entier */
      *param_int = stringToInt(ligne);
    
    else if (cpt==6 || cpt==7){ 
      stringToString(*(param_char[cpt_char]),ligne);
      cpt_char++;
    }

    else{
      /* on doit convertir la chaine en double */      
      *(param_double[cpt_double]) = stringToDouble(ligne);      
      cpt_double++;
    }

  }

  fclose(fic_in_impsig);

  free(ligne);

}

void loadParamOptim(char *name_in_optim, double *pt_gradtol, double *pt_steptol, int *pt_verbosity, int *pt_saveSuccessiveXinFile, int *pt_maxCounter, double *pt_lambda){

  char *ligne;
  char *ok;
  int nbr;
  int cpt;
  double *param_double[3];
  int *param_int[3];

  int cpt_double;
  int cpt_int;

  int i=0;

  FILE *fic_in_optim;
  fic_in_optim = fopen(name_in_optim,"r");

  param_double[0] = pt_gradtol;
  param_double[1] = pt_steptol;
  param_double[2] = pt_lambda;

  param_int[0] = pt_verbosity;
  param_int[1] = pt_saveSuccessiveXinFile;
  param_int[2] = pt_maxCounter;

  cpt_double = 0;
  cpt_int = 0;

  nbr = 100;
  
  ligne = (char *) malloc(nbr*sizeof(char));

  for (cpt=1;cpt<=6;cpt++){

    ok = fgets(ligne,nbr,fic_in_optim);
    if (ok == NULL){
      printf("Pb de lecture dans le fichier d'entree\n");
      exit;
    }
    
    while (ligne[0] == '#'){
      ok = fgets(ligne,nbr,fic_in_optim);
      if (ok == NULL){
	printf("Pb de lecture dans le fichier d'entree\n");
	exit;
      }
    }
    
    /* dans ligne : la valeur numero cpt+1 sous forme de chaine de carac */
    if (cpt==3 || cpt==4 || cpt==5){
      /* on doit convertir la chaine en entier */
      *(param_int[cpt_int]) = stringToInt(ligne);
      cpt_int++;
    }
    else{
      /* on doit convertir la chaine en double */      
      *(param_double[cpt_double]) = stringToDouble(ligne);      
      cpt_double++;
    }

  }

  fclose(fic_in_optim);

  free(ligne);

}

void loadParamOptim2(char *name_in_optim, double *pt_pgtol, double *pt_factr, int *pt_iprint, int *pt_maxCounter, double *pt_sigma_min, double *pt_sigma_max, double *pt_lambda){

  char *ligne;
  char *ok;
  int nbr;
  int cpt;
  double *param_double[5];
  int *param_int[2];

  int cpt_double;
  int cpt_int;

  int i=0;

  FILE *fic_in_optim;
  fic_in_optim = fopen(name_in_optim,"r");

  param_double[0] = pt_pgtol;
  param_double[1] = pt_factr;
  param_double[2] = pt_sigma_min;
  param_double[3] = pt_sigma_max;
  param_double[4] = pt_lambda;


  param_int[0] = pt_iprint;
  param_int[1] = pt_maxCounter;

  cpt_double = 0;
  cpt_int = 0;

  nbr = 100;
  
  ligne = (char *) malloc(nbr*sizeof(char));

  for (cpt=1;cpt<=7;cpt++){

    ok = fgets(ligne,nbr,fic_in_optim);
    if (ok == NULL){
      printf("Pb de lecture dans le fichier d'entree\n");
      exit;
    }
    
    while (ligne[0] == '#'){
      ok = fgets(ligne,nbr,fic_in_optim);
      if (ok == NULL){
	printf("Pb de lecture dans le fichier d'entree\n");
	exit;
      }
    }
    
    /* dans ligne : la valeur numero cpt+1 sous forme de chaine de carac */
    if (cpt==3 || cpt==4){
      /* on doit convertir la chaine en entier */
      *(param_int[cpt_int]) = stringToInt(ligne);
      cpt_int++;
    }
    else{
      /* on doit convertir la chaine en double */      
      *(param_double[cpt_double]) = stringToDouble(ligne);      
      cpt_double++;
    }

  }

  fclose(fic_in_optim);

  free(ligne);

}

void stringToString(char *strcopy, char *string){

  int size=0;
  int i=0;

  while (string[i] != '\n'){
    size ++;
    i++;
  }
  size --;

  for (i=0;i<=size;i++)
    strcopy[i] = string[i];
  strcopy[size+1]='\0';
  
}

int stringToInt(char *string){

  /*string se termine par \n\0 */
  int res=0;
  int i=0;
  int size=0;

  while (string[i] != '\n'){
    size ++;
    i++;
  }
  size --;

  for (i=0;i<=size;i++)
    res = res + (string[i]-'0') * pow(10,size-i);

  //printf("---> %s + taille %d \n",string,size);

  return res;

}

double stringToDouble(char *string){

  int index_point=0; /* indice de position du point, = 1000 si pas de point */
  int size=0;
  char *string2;
  double res=0;
  int i;

  int nbr = 100;

  string2 = (char *) malloc(nbr*sizeof(char));

  if (string[0] == '-'){
    i=0;
    while (string[i+1] != '\n'){
      string2[i] = string[i+1];
      i++;
    }
    string2[i] = '\n';
  }
  else{
    i=0;
    while (string[i] != '\n'){
      string2[i] = string[i];
      i++;
    }
    string2[i] = '\n';
  }

  /* dans string2 : la chaine a traiter (representant un double positif) */
  i = 0;
  while (string2[i] != '\n'){
    size ++;
    i++;
  }
  size--;

  while (string2[index_point] != '\n' && string2[index_point] != '.' )
    index_point ++;

  if (string2[index_point] == '\n'){
    string2[size+1] = '.';
    string2[size+2] = '\n';
    size++;
    /* index_point  == size */
  }

  
  for (i=0;i<index_point;i++)
    res = res + (string2[i]-'0')*pow(10,index_point-1-i);

  for (i=index_point+1;i<=size;i++)
    res = res + (string2[i]-'0')*pow(10,index_point-i);

  free(string2);

  if (string[0] == '-')
    return -res;
  else
    return res;

}

void savePriceOrSigma(char *name_out_visu, double **priceOrSigma_visuGrid, int Nprice_visu, int Mprice_visu, double *Kprice_visu, double *Tprice_visu){

  int i,j;
  FILE *fic_sortie;

  fic_sortie = fopen(name_out_visu,"w");

  fprintf(fic_sortie,"%d\n",Nprice_visu);
  fprintf(fic_sortie,"%d\n",Mprice_visu);
  for (i=0;i<=Nprice_visu;i++)
      fprintf(fic_sortie,"%lf\n",Kprice_visu[i]);
  for (j=0;j<=Mprice_visu;j++)
      fprintf(fic_sortie,"%lf\n",Tprice_visu[j]);
  for (i=0;i<=Nprice_visu;i++)
    for (j=0;j<=Mprice_visu;j++)
	fprintf(fic_sortie,"%lf\n",priceOrSigma_visuGrid[i][j]);

  fclose(fic_sortie);

}

void saveSigmaddl(char *name_out_sigma_ddl, double *sigma_ddl, int n, int m, double *y_coarseGrid, double *T_coarseGrid){

  int i,nbParam;
  FILE *fic_out_sigma_ddl;

  nbParam = (n+3)*(m+3); // i.e (n+1)*(m+1)+2*(m+1)+2*(n+1)+4 

  fic_out_sigma_ddl = fopen(name_out_sigma_ddl,"w");

  fprintf(fic_out_sigma_ddl,"%d\n",n);
  fprintf(fic_out_sigma_ddl,"%d\n",m);
  for (i=0;i<=n;i++)
      fprintf(fic_out_sigma_ddl,"%lf\n",y_coarseGrid[i]);
  for (i=0;i<=m;i++)
      fprintf(fic_out_sigma_ddl,"%lf\n",T_coarseGrid[i]);
  for (i=0;i<nbParam;i++)
      fprintf(fic_out_sigma_ddl,"%lf\n",sigma_ddl[i]);

  fclose(fic_out_sigma_ddl);

}
