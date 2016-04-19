#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "tree.h"
#include "rollback.h"
#include "inout.h"

#include "optim-code/QuasiNewton.h"

Trinomial_tree arbre ;
Tri_tree_prob_space V ;
Tri_tree_prob_space GV ;
Vecteur indice , diffMatu ;
Vect_option SOption ;
Parametre param ;
double *lamda;
double *source ;
double *Gsource;
double alpha;
double priorInput;


int calib_ave_scilab ()
{
	int i , j = 0 ;
	int N, max_counter,verbosity,saveSuccessiveXinFile;
	double grad_tol , step_tol ,lambda0;
	double S_0,r,q,sigma_0,sigma_min,sigma_max,sigma_bar;

	double prixOptim,erreurPrixOptim=0.0 , impVol, maxVol , minVol , tempVol , erreur=0.0 , dt,t;
  	double T_max;

  	double *grad, *sol ;

	char *name_in_data;
        char *name_out_Vol_Loc;
	char *name_in_calib="calib_Avellaneda.in";
	char *name_in_pricer="pricer_Avellaneda.in";


	Vect_option option ;


	FILE *ArbreFile;
	FILE *Pricer_in;
        FILE *VolLoc;
	alpha=1;


	loadParamCalib(name_in_calib, &S_0, &r,&q, &N , &sigma_0, &sigma_min, &sigma_max, &sigma_bar,
&grad_tol,&step_tol,&verbosity,&saveSuccessiveXinFile, &max_counter,&lambda0,&alpha,&name_in_data,&name_out_Vol_Loc);

	param.S0 = S_0;
  	param.N = N ;
	param.dividende = q;
  	param.r = r;
	param.mu = param.r - param.dividende;
	param.prior = sigma_0;
	param.prior_max = sigma_max ;
  	param.prior_min = sigma_min ;
	param.sigmabar =sigma_bar;


	printf("alpha %f lamda %f prior %f\n",alpha,lambda0,sigma_0);

	//chargement des donnees dans la structure option
	option=load_options( name_in_data );

	// tri des donnees par maturite croissante dans la structure SOption (Sorted Option)
	SOption = Sort(option);

	//renvoie uniquement les DIFFERENTES maturites dans l'ordre croissant.
  	diffMatu = getDiffMaturity(option) ;

	// IL FAUT TROUVER LE PLUS PETIT dt  POUR QUE LES MATURITES COINCIDENT AVEC L'ARBRE TRINOMIAL
  	//ON NE DOIT PAS CHANGER LE PAS DE TEMPS EN COURS DE ROUTE SINON L'ARBRE N'EST PLUS RECOMBINANT
  	dt=getdt(diffMatu,3);
	T_max = SOption.maturite[SOption.size-1];
	param.dt = T_max/N;

  	param.dt = diffMatu.val[diffMatu.size-1]/(double)param.N ;

	impVol = 0.0;
	minVol=10;
	maxVol=0.0;

	for(i=0;i<SOption.size; i++)
	{
      tempVol = ImpliedVol(param.S0, param.r, param.dividende, &SOption.payoff[i] , 0 , SOption.strik[i], SOption.maturite[i], SOption.prix[i],0.000001);
      impVol += tempVol ;
      if (tempVol<minVol) minVol=tempVol ;
      if (tempVol>maxVol) maxVol=tempVol ;
      //prixBS = BS(&SOption.payoff[i],param.S0, tempVol , param.r , SOption.strik[i] , SOption.maturite[i],param.dividende);
      printf("T %0.3f K %3f vol implicite %f prixfichier %f Flag %c\n",SOption.maturite[i],SOption.strik[i] ,tempVol,SOption.prix[i],SOption.payoff[i]);

	}
	param.prior = impVol/SOption.size ;
	printf("prior %f MaxImpVol %f MinImpVol %f\n",param.prior, maxVol,minVol);

	printf("prior %f priorMax %f sigmabar %f priorMin %f\n",param.prior, param.prior_max,param.sigmabar,param.prior_min);

	indice.size = diffMatu.size ;
	indice.val = (double *) malloc(diffMatu.size*sizeof(double));
	arbre = getPriorTree( diffMatu, &indice , param);

	lamda =(double *)malloc(SOption.size*sizeof(double));
	for(i=0 ; i< SOption.size ; i++){lamda[i]=lambda0;}

	V = define_3_tree_prob_space( arbre.tdim , 0.0 );
	GV = define_3_tree_prob_space( arbre.tdim , 0.0 );

	grad =(double *)malloc(SOption.size*sizeof(double)) ;

	sol =(double *)malloc(SOption.size*sizeof(double)) ;

	source = (double *)calloc((arbre.tdim+1)*(arbre.tdim+1),sizeof(double));
	Gsource = (double *)calloc((arbre.tdim+1)*(arbre.tdim+1),sizeof(double));

	// debut de QuasiNewton //

    printf("avant QuasiNewton\n");

    QuasiNewton( SOption.size ,
	                lamda,
	                costFunction,
	               	gradCostFunction,
	                sol,
	                grad_tol,
	                step_tol,
	                max_counter,
	                verbosity,
	                saveSuccessiveXinFile);

	printf("Lagrangien : lamda optimal \n");
	for( i = 0 ; i<SOption.size ;i++){printf("sol %f \n",sol[i]);}
	printf("fin sol\n");


	printf("K  T  C  Cmodel  ErrRel  VolImpModel\n");


	for( i = 0 ; i<SOption.size ;i++)
	{

		prixOptim=optimPrice(SOption.strik[i],SOption.maturite[i],param,&SOption.payoff[i]);
		tempVol = ImpliedVol(param.S0, param.r, param.dividende, &SOption.payoff[i] , 0 , SOption.strik[i], SOption.maturite[i], prixOptim,0.00000001);
		erreur =fabs(prixOptim-SOption.prix[i])/SOption.prix[i];
		//printf("prix initial %f prix %f erreur %f impVol %f \n",SOption.prix[i] ,prixOptim, erreur, tempVol );
		printf("%f %f %f %f %f %f\n",SOption.strik[i],SOption.maturite[i],SOption.prix[i],prixOptim,erreur,tempVol);
		erreurPrixOptim +=erreur;
	}
	printf("erreurPrixOptim %10.10e \n",erreurPrixOptim);

	//printf("gradient optimal\n");
	//gradCostFunction ( sol , grad );

	for( i = 0 ; i<SOption.size ;i++){printf("grad %f \n",grad[i]);}
	ArbreFile=fopen("arbre_Avellaneda.out","w");
        VolLoc=fopen(name_out_Vol_Loc,"w");
	fprintf(ArbreFile,"#############################################################################\n") ;
	fprintf(ArbreFile,"#####################################S_O#####################################\n") ;
	fprintf(ArbreFile,"%f\n",param.S0) ;
	fprintf(ArbreFile,"#####################################r:taux constant#########################\n") ;
	fprintf(ArbreFile,"%f\n",param.r) ;
	fprintf(ArbreFile,"#####################################q:dividende#############################\n") ;
	fprintf(ArbreFile,"%f\n",param.dividende) ;
	fprintf(ArbreFile,"#####################################N: discretisation en temps##############\n") ;
	fprintf(ArbreFile,"%d\n",(int)param.N) ;
	fprintf(ArbreFile,"#####################################sigma_0: prior utilise dans l'arbre#####\n") ;
	fprintf(ArbreFile,"%f\n",param.prior) ;
	fprintf(ArbreFile,"###############################sigma_min: sigma_min utilise dans l'arbre#####\n") ;
	fprintf(ArbreFile,"%f\n",param.prior_min) ;
	fprintf(ArbreFile,"###############################sigma_max: sigma_max utilise dans l'arbre#####\n") ;
	fprintf(ArbreFile,"%f\n",param.prior_max) ;
	fprintf(ArbreFile,"###############################sigma_bar: sigma_bar utilise dans l'arbre#####\n") ;
	fprintf(ArbreFile,"%f\n",param.sigmabar) ;
	fprintf(ArbreFile,"###############################T: Temps final de l'arbre#####################\n") ;
	fprintf(ArbreFile,"%f\n",diffMatu.val[diffMatu.size-1]) ;
        fprintf(ArbreFile,"#####################Input File for the calibrated Trinomial Tree############\n") ;
	fprintf(ArbreFile,"%s\n",name_in_data) ;
	fprintf(ArbreFile,"####################################ARBRE TRINOMIAL##########################\n") ;
	fprintf(ArbreFile,"#####################################S dt proba##############################\n") ;
	fprintf(ArbreFile,"#############################################################################\n") ;

        fprintf(VolLoc,"#####################################S_O#####################################\n") ;
	fprintf(VolLoc,"%f\n",param.S0) ;
	fprintf(VolLoc,"#####################################r:taux constant#########################\n") ;
	fprintf(VolLoc,"%f\n",param.r) ;
	fprintf(VolLoc,"#####################################q:dividende#############################\n") ;
	fprintf(VolLoc,"%f\n",param.dividende) ;
	fprintf(VolLoc,"#####################################N: discretisation en temps##############\n") ;
	fprintf(VolLoc,"%d\n",(int)param.N) ;
	fprintf(VolLoc,"#####################################sigma_0: prior utilise dans l'arbre#####\n") ;
	fprintf(VolLoc,"%f\n",param.prior) ;
	fprintf(VolLoc,"###############################sigma_min: sigma_min utilise dans l'arbre#####\n") ;
	fprintf(VolLoc,"%f\n",param.prior_min) ;
	fprintf(VolLoc,"###############################sigma_max: sigma_max utilise dans l'arbre#####\n") ;
	fprintf(VolLoc,"%f\n",param.prior_max) ;
	fprintf(VolLoc,"###############################sigma_bar: sigma_bar utilise dans l'arbre#####\n") ;
	fprintf(VolLoc,"%f\n",param.sigmabar) ;
	fprintf(VolLoc,"###############################T: Temps final de l'arbre#####################\n") ;
	fprintf(VolLoc,"%f\n",diffMatu.val[diffMatu.size-1]) ;
        fprintf(VolLoc,"#####################Input File for the calibrated Trinomial Tree############\n") ;
	fprintf(VolLoc,"%s\n",name_in_data) ;
	fprintf(VolLoc,"####################################ARBRE TRINOMIAL##########################\n") ;
	fprintf(VolLoc,"#####################################S t VolLoc##############################\n") ;
	t=0;
	for ( i = 0 ; i < arbre.tdim + 1 ; i++)
 	{
                for ( j = i*i ; j < (i+1)*(i+1) ; j++)
                {
                        fprintf(VolLoc,"%f %f %f \n",arbre.val[j],t ,param.sigmabar*sqrt(arbre.proba[j])) ;
                        fprintf(ArbreFile,"%f %f %f \n",arbre.val[j], arbre.pasTemps[j],arbre.proba[j]) ;
                }
                t+=param.dt;
  	}

	fclose(ArbreFile);
        fclose(VolLoc);
	Pricer_in=fopen(name_in_pricer,"w");
	fprintf(Pricer_in,"#############################################################################\n") ;
	fprintf(Pricer_in,"##############Used Parameters To calibrate the Tree #########################\n") ;
	fprintf(Pricer_in,"########S_O= %.3f r= %.3f dividende q = %.3f N=%.0f #####################\n", param.S0,param.r,param.dividende,param.N);
	fprintf(Pricer_in,"#####prior sigma_0=%.3f sigma_min=%.3f sigma_min =%.3f sigmabar=%.3f ####\n", param.prior,param.prior_min,param.prior_max,param.sigmabar);
	fprintf(Pricer_in,"############################################################\n") ;
	fprintf(Pricer_in,"###################### VARIABLES K, T, AND OPTIONTYPE #######################\n") ;
	fprintf(Pricer_in,"#############################################################################\n") ;
	fprintf(Pricer_in,"# K : Strike of the option to price\n") ;
	fprintf(Pricer_in,"%f\n",SOption.strik[SOption.size-1]) ;
	fprintf(Pricer_in,"# T : Maturity of the option to price must be less then the \n") ;
	fprintf(Pricer_in,"# terminal date on the calibrated trinomial Tree =%.3f\n",SOption.maturite[SOption.size-1]) ;
	fprintf(Pricer_in,"%f\n",T_max) ;
	fprintf(Pricer_in,"# optionType : type of the option (C for call, P for put)\n") ;
	fprintf(Pricer_in,"%c\n",SOption.payoff[SOption.size-1]) ;
	fprintf(Pricer_in,"#############################################################################\n") ;
        fprintf(Pricer_in,"############file containing the option to price##############################\n") ;
        fprintf(Pricer_in,"##############the data must be stored in this order #########################\n");
        fprintf(Pricer_in,"#################### K T flag(P as put,Cas Call) ############################\n");
        fprintf(Pricer_in,"\n") ;
        fprintf(Pricer_in,"############output file containing the price of these options################\n") ;
        fprintf(Pricer_in,"\n") ;
        fprintf(Pricer_in,"#############################################################################\n") ;
        fclose(Pricer_in);

        free(sol) ;
	free(grad) ;

        free(name_in_data);
        free(name_out_Vol_Loc);

	free_Trinomial_tree( arbre ) ;
  	free_Vect_option(option) ;
  	free_Vect_option(SOption) ;
  	free(lamda) ;
	free(source);
  	free(Gsource);
  	free_Tri_tree_prob_space(V);
	free_Tri_tree_prob_space(GV);
	free_Vecteur(indice) ;
	free_Vecteur(diffMatu) ;

	return 0;
}
