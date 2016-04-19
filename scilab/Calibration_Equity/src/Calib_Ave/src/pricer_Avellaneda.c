/***************************************************************************
                          pricer_Avellaneda.c  -  description
                             -------------------
    begin                : Thu Jan 15 2004
    copyright            : (C) 2004 by messaoud
    email                : marouen.messaoud@inria.fr
 ***************************************************************************/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tree.h"
#include "inout.h"
#include "pricer_Avellaneda.h"

Trinomial_tree arbre ;

Parametre paramPricer;

void loadTree(char *name_in_tree , double *pt_S_0, double *pt_r, double *pt_q , int *pt_N ,
	double *pt_sigma_0, double *pt_sigma_min, double *pt_sigma_max, double *pt_sigma_bar, double *pt_T_max,char **pt_name_in_data)
{
	char *ligne;
    char *ok;
    double s,p,dt;
    int cpt;
    double *param_double[8];
    int *param_int[1];
    char **param_char[1];

    int cpt_double;
    int cpt_int;
    int cpt_char;
    int nbr;
    int i=0;

    FILE *fic_in_tree;
    fic_in_tree = fopen(name_in_tree,"r");

    param_double[0] = pt_S_0;
    param_double[1] = pt_r;
    param_double[2] = pt_q;
    param_double[3] = pt_sigma_0;
    param_double[4] = pt_sigma_min;
    param_double[5] = pt_sigma_max;
    param_double[6] = pt_sigma_bar;
    param_double[7] = pt_T_max;

    param_char[0] = pt_name_in_data;
    nbr = 100;

    param_int[0] = pt_N;

    cpt_double = 0;
    cpt_int = 0;
    cpt_char = 0;

    nbr = 100;
    //1 comme deux parametre char**
    for (i=0;i<1;i++)
        *(param_char[i]) = (char *) malloc(nbr*sizeof(char));

    ligne = (char *) malloc(nbr*sizeof(char));
	// 10 parametres
	for (cpt=1;cpt<=10;cpt++)
  	{

    	ok = fgets(ligne,nbr,fic_in_tree);
	    if (ok == NULL)
	    {
    		printf("Pb de lecture dans le fichier d'entree\n");
	      	exit(-1);
    	}

    	while (ligne[0] == '#')
	{
			ok = fgets(ligne,nbr,fic_in_tree);
      		if (ok == NULL)
	      	{
				printf("Pb de lecture dans le fichier d'entree\n");
				exit(-1);
	        }
    	}
        if (cpt==4){
	      	/* on doit convertir la chaine en entier */
    	  	*(param_int[cpt_int]) = stringToInt(ligne);
      		cpt_int++;
    	}
        else if (cpt==10){
                stringToString(*(param_char[cpt_char]),ligne);
                cpt_char++;
      }
		else
		{
	    	/* dans ligne : la valeur numero cpt+1 sous forme de chaine de carac */
	    	/* on doit convertir la chaine en double */
		    *(param_double[cpt_double]) = stringToDouble(ligne);
    		cpt_double++;
    	}
	}

	arbre=define_3_tree((int) *pt_N , 0.0 );

 	for ( i = 0 ; i < ( (int)*pt_N+1)*((int)*pt_N+1)  ; i++)
	{
		ok = fgets(ligne,nbr,fic_in_tree);
		if (ok == NULL)
	    {
    		printf("Pb de lecture dans le fichier d'entree\n");
	      	exit(-1);
    	}
		while (ligne[0] == '#')
	    {
			ok = fgets(ligne,nbr,fic_in_tree);
      		if (ok == NULL)
	      	{
				printf("Pb de lecture dans le fichier d'entree\n");
				exit(-1);
			}
    	}

		sscanf(ligne,"%lf %lf %lf\n", &s,&dt,&p) ;

	    arbre.val[i] = s ;
	    arbre.pasTemps[i] = dt ;
	    arbre.proba[i] = p ;

	}
        fclose(fic_in_tree);
	free(ligne);
}

double optimPrice(double K, double T, Parametre param , char *OptionType)
{
		Tri_tree_prob_space Vprix ;
		int i,j = 0 ,n, dim = 1;
		double s , temp , Pm , Pu , Pd , temp0 , t , p , dt ;
		char message[80];
		Vprix = define_3_tree_prob_space( arbre.tdim , 0.0 );

		sprintf(message,"%f",T/param.dt);
		sscanf(message,"%i",&n);

		t=T;
		for ( j = n*n ; j < (n+1)*(n+1) ; j++)
	    {
    		s = arbre.val[j];
      		Vprix.val[0][j] = Gi(s,K,OptionType);
		}

	    for ( i = n-1 ; i >= 0 ; i--)
    	{
     		dt = arbre.pasTemps[i*i] ;
     		t-= dt ;

      		for ( j = i*i ; j < (i+1)*(i+1) ; j++)
      		{
				s = arbre.val[j] ;
		  		temp0 = param.sigmabar*sqrt(dt)/2.0 ;

		  		p=arbre.proba[j];

			  	//p = 1 ;
			  	Pu = 0.5*p*(1-temp0) +0.5*param.mu*sqrt(dt)/param.sigmabar ;
        		Pm = 1-p ;
        		Pd = 0.5*p*(1+ temp0) -0.5*param.mu*sqrt(dt)/param.sigmabar ;
			  	Vprix.val[0][j] = Pu* Vprix.val[0][2*i+j+3]+Pm*Vprix.val[0][2*i+j+2]+Pd*Vprix.val[0][2*i+j+1] ;
			  	Vprix.val[0][j] = Vprix.val[0][j];
			   //printf("somme des probas %f\n",Pu+Pm+Pd);
			 }
		}
	    temp = Vprix.val[0][0]*exp(-param.r*T) ;

		free_Tri_tree_prob_space(Vprix) ;

		return(temp) ;
}

int main()
{
	int N;

	double S_0,r,q,sigma_0,sigma_min,sigma_max,sigma_bar,T_max;
	double K,T;
	double prixOptim;
    double Ki,Ti ;
	//char payoff ;
	//char *flag;
	int i = 0 , j = 0 ;
	char ch_line[200];
	FILE *fi ;
    FILE *fout ;

  	char *option_type;
    char *name_in_grid_price_in;
    char *name_in_grid_price_out;
    char *name_in_data;
	char *name_in_pricer="pricer_Avellaneda.in";
	char *name_in_tree="arbre_Avellaneda.out";

	loadTree(name_in_tree,&S_0, &r,&q ,&N ,&sigma_0,&sigma_min,&sigma_max,&sigma_bar,&T_max,&name_in_data);

	paramPricer.S0 = S_0;
  	paramPricer.N = N ;
	paramPricer.dividende = q;
  	paramPricer.r = r;
	paramPricer.mu = paramPricer.r - paramPricer.dividende;
	paramPricer.prior = sigma_0;
	paramPricer.prior_max = sigma_max ;
  	paramPricer.prior_min = sigma_min ;
	paramPricer.sigmabar =sigma_bar;
	paramPricer.dt = T_max/N;

	loadParamPricer(name_in_pricer,&K, &T, &option_type,&name_in_grid_price_in,&name_in_grid_price_out);

        if(K>=0 && T>=0){
                if( T>T_max){printf("the maturity choosed is greater then the final time of the calibrated tree %f\n",T_max);return 0;}
                else{
	                prixOptim=optimPrice(K,T,paramPricer,option_type);
                        printf("K %f T %f Price %f optionType %s\n",K, T, prixOptim,option_type) ;
                }
        }else{
                printf("problem with the parametre\n");
                return 0;
        }

        if ( *name_in_grid_price_in == 0 && *name_in_grid_price_out != 0)
		{
			printf("YOU MUST CHOOSE AN INPUT FILE IN pricer_Avellaneda.in Or remove the output file name\n");
			return 0;
		}
		else if ((fi=fopen(name_in_grid_price_in,"r")) != NULL)
		{
		    while( fgets(ch_line,200,fi)!= NULL ) j++;
   	        fclose(fi);
	        printf("Nombre d'option dans le fichier %s = %d\n",name_in_grid_price_in,j);
            //fin get number of line
            fout=fopen(name_in_grid_price_out,"w");

            fprintf(fout,"Prices Computed with the last output of calib_Avellaneda\n(file %s) from the data file %s\n",name_in_tree,name_in_data) ;
            fprintf(fout,"K T Price OptionType\n") ;

            while((fi=fopen(name_in_grid_price_in,"r"))==NULL)
		    printf("impossible d'ouvrir le fichier %s \n",name_in_grid_price_in);


	        while(fscanf(fi,"%lf %lf %s\n", &Ki,&Ti,option_type)!=EOF)
	        {
	                prixOptim=optimPrice(Ki,Ti,paramPricer,option_type);

                    fprintf(fout,"%f %f %f %s\n",Ki, Ti, prixOptim,option_type) ;
             }

	        fclose(fi);
            fclose(fout);
        }
        //printf("price from calibrated Tree %f %s %s \n",prixOptim,name_in_grid_price_in,name_in_grid_price_in );
        free(option_type);
        free(name_in_grid_price_in);
        free(name_in_grid_price_out);
        free(name_in_data);

	return 0;
}
