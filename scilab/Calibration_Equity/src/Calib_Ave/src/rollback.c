#include "tree.h"
#include "rollback.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern Trinomial_tree arbre ;
extern Tri_tree_prob_space V ;
extern Tri_tree_prob_space GV ;
extern Vecteur indice , diffMatu ;
extern Vect_option SOption ;
extern Parametre param ;
extern double *lamda;
extern double alpha;
extern double *source ;
extern double *Gsource;

double PHI( double x , double prior , double minp , double maxp )
{
  double res = 0.0 ;
	double low = 2*alpha*log(minp/prior);
	double hy = 2*alpha*log(maxp/prior);
	
	if ( x < hy && x > low )
  {
    res = alpha*prior*prior*(exp(x/alpha)-1) ;
  }
  else if ( x <= low )
  {
    res = minp*minp*x -minp*minp*low+alpha*minp*minp-alpha*prior*prior;
  }
  else
  {
    res = maxp*maxp*x-maxp*maxp*hy+alpha*maxp*maxp-alpha*prior*prior;
  }
	/*
  double low = 2*alpha*(minp*minp - prior*prior) ;
  double hy = 2*alpha*(maxp*maxp - prior*prior) ;

  if ( x < hy && x > low )
  {
    res = 0.25*x*x/alpha+prior*prior*x ;
  }
  else if ( x <= low )
  {
    res = minp*minp*x - 0.25*(low*low)/alpha ;
  }
  else
  {
    res = maxp*maxp*x - 0.25*(hy*hy)/alpha ;
  }*/

  return (res) ;  
}

double GPHI( double x , double prior , double minp , double maxp )
{
  double res = 0.0 ;
	double low = 2*alpha*log(minp/prior);
	double hy = 2*alpha*log(maxp/prior);
	
	if ( x < hy && x > low )
  {
    res = prior*prior*exp(x/alpha) ;
  }
  else if ( x <= low )
  {
    res = minp*minp;
  }
  else
  {
    res = maxp*maxp;
  }

  /*double low = 2*alpha*(minp*minp - prior*prior) ;
  double hy = 2*alpha*(maxp*maxp - prior*prior) ;

  if ( x < hy && x > low )
  {
    res = 0.5*x/alpha + prior*prior ;
  }
  else if ( x <= low )
  {
    res = minp*minp ;
  }
  else
  {
    res = maxp*maxp ;
  }*/

  return (res) ;
}

Trinomial_tree getPriorTree(/*vecteur des differentes maturites ordonnees */ Vecteur diffMaturite,Vecteur *indice , Parametre param )
{
    int  k=0, i , j , N ;
    double Ti , t , S0=param.S0 , sigmabar = param.sigmabar ;
    double dt = param.dt, sdt , temp ;
    Vecteur vect ;

    Trinomial_tree Tree;
    N = (int) param.N ;

    vect.size = N+1 ;
    vect.val = (double *) malloc(vect.size*sizeof(double));

    t = 0 ;
    j = 0 ;
    for ( i = 0 ; i < diffMaturite.size ; i++)
    {
      Ti = diffMaturite.val[i] ;
      while( t  <  Ti -0.00001)
      {
        vect.val[k] = t ;
        k++;
        t+= dt ;
      }
      if ( fabs(t -Ti)< 0.00001 )
      {
        vect.val[k] = Ti ;
        (*indice).val[j] = (double)k ;
        j++ ;
        k++ ;
        t+= dt ;
      }

    }

	Tree=define_3_tree( N , 0.0 );

    Tree.val[0] = param.S0 ;
    Tree.pasTemps[0] = dt ;


    for ( i = 1 ; i < N+1 ; i++)
    {
      sdt = sqrt((double)dt) ;
      temp = exp(sigmabar*sdt) ;
      Tree.val[i*i] = Tree.val[(i-1)*(i-1)]/temp;
      Tree.pasTemps[i*i] = dt ;

      for ( j = i*i +1 ; j < (i+1)*(i+1) ; j++)
      {
        Tree.pasTemps[j] = dt ;
        Tree.val[j] = Tree.val[j-1] * temp ;
      }

    }
	free_Vecteur(vect) ;
	return Tree ;


}


void define_G_source(double *Gsource , double Ki , char *payoff  , Trinomial_tree arbre , int rank  )
{

    int j ;
    double  S;

   	for( j = rank*rank ; j < (rank+1)*(rank+1) ; j++ )
  	{
		S = arbre.val[j] ;

   		Gsource[j]+= Gi( S , Ki , payoff);
   	}

}

/* define_source retourne le vecteur source sous la structure arbre a partir
de la collection des options et de l'arbre discretise de l'espace*/
int define_source( double *source , Vect_option SOption , Trinomial_tree arbre ,Vecteur indice , double *lam , double r )
{
    int dim = SOption.size ;

    int n , m , i , j ;
    double K , C , Lam ,Ti , S;

    m = 0 ;

    for( i = 0 ; i < dim ;  i++ )
    {
    	Ti = SOption.maturite[i] ;
    	n = (int)indice.val[m] ;
    	
    	/*l`indice de ces options de mat Ti est dans indice
        le vecteur option a ete ordonne selon les Ti croissant*/
		K = SOption.strik[i] ;
 		C = SOption.prix[i] ;
        Lam = lam[i] ;

       	for( j = n*n ; j < (n+1)*(n+1) ; j++ )
        {
           	S = arbre.val[j] ;
            source[j]+= Lam * ( Gi( S , K ,&SOption.payoff[i]) ) ;
   	   	}
        m=m+(Ti!=SOption.maturite[i+1]) ;

    }
		
	return 1 ;

}


double costFunction ( double *lam )
{

    double p , Pu , Pm , Pd ;
    double r , mu , prior , sigmabar , maxprior , minprior , dt=0.0 , S , t ;
	double gama ; // derive seconde de V %S  (Vss)
	double temp = 0.0 , temp0 =0.0 , temp1 = 0.0;
	int n , i , j , tdim ;
	
	tdim = arbre.tdim ;
	r = param.r ;
	mu = param.mu ;
	prior = param.prior ;
	sigmabar = param.sigmabar ;
	minprior = param.prior_min ;
	maxprior = param.prior_max ;
	n = V.tdim ;// le nombre de date dans V, (n+1)^2=taille de V
	

  	for( i =0 ; i < (tdim+1)*(tdim+1) ; i++)
  	{
  		source[i] = 0.0;
  	}
	
  	define_source( source , SOption , arbre , indice , lam , r ) ;
		
	t = SOption.maturite[SOption.size-1] ; // t = T
		
  	i = n ;

	for ( j = i*i ; j < (i+1)*(i+1) ; j++)
  	{
  		V.val[0][j] = source[j] ;
	}

 	for ( i = n-1 ; i >= 0 ; i--)
  	{
 		dt = arbre.pasTemps[i*i] ;
     	t-= dt ;
			
      	for ( j = i*i ; j < (i+1)*(i+1) ; j++)
      	{
				
			// calcul du parametre p(i,j) , Gp(i,j) a partir de phi , phi' et de V a la date i+1 et au niveau (j,j+1,j+2)
        	//p(i,j)=(phi(exp(-0.5rt S^2 Vss)/(sigmabar^2 exp(-rt) 0.5 S^2 Vss)
        	// gama = 0.5 S^2 Vss
        	S = arbre.val[j] ;
			temp0 = sigmabar*sqrt(dt)/2.0 ;
	        gama = (( 1 - temp0 )*V.val[0][2*i+j+3]+( 1.0 + temp0 )*V.val[0][2*i+j+1]- 2.0*V.val[0][2*i+j+2] )/(4*temp0*temp0);
    	    temp = exp(-r*t)*gama ;
				
        	if (fabs(temp)>1e-8) p = PHI( temp , prior , minprior , maxprior )/ ( temp*sigmabar*sigmabar ) ;
        	else p = prior*prior/(sigmabar*sigmabar) ;
        	
			V.val[1][j] = p ;

			arbre.proba[j]=GPHI(temp,prior,minprior,maxprior)/(sigmabar*sigmabar);
    		
        	//calcul des proba Pu(i,j) , Pm(i,j), Pd(i,j) et de  GPu(i,j) , GPm(i,j), GPd(i,j)
        	//Pu=0.5p(1-0.5sigmabar sqrt(dt) ) +0.5 mu sqrt(dt)/sigmabar
        	//Pm=1-p
        	//Pd= 0.5  p(1+0.5sigmabar sqrt(dt) ) -0.5 mu sqrt(dt)/sigmabar
        	Pu = 0.5*p*(1-temp0) +0.5*mu*sqrt(dt)/sigmabar ;
	        Pm=1-p;
    	    Pd= 0.5*p*(1+ temp0) -0.5*mu*sqrt(dt)/sigmabar ;
				
			if(Pu<=0 || Pu>=1  || Pm>=1 || Pm<=0 || Pd>=1|| Pd<=0  )
			{ printf("\n\n\nPROBA NEGATIVE dans V PU %f Pm %f Pd %f \n\n",Pu,Pm,Pd);}
    	    //calcul de V(i,j)=Pu(i,j)V(i+1,j+2)+Pm(i,j)V(i+1,j+1)+Pd(i,j)V(i+1,j)+source(i,j)
        	V.val[0][j] = Pu* V.val[0][2*i+j+3]+Pm*V.val[0][2*i+j+2]+Pd*V.val[0][2*i+j+1] ;
	        V.val[0][j] = V.val[0][j] * exp(-r*dt) + source[j];
		}
	}
	
	for (i = 0 ; i < SOption.size ; i++){V.val[0][0]-=lam[i]*SOption.prix[i] ;}

	return V.val[0][0];

}


void gradCostFunction ( double *lam , double *grad )
{

	int k , rank , m = 0;
	int n , i , j , tdim ;

	double Gp, GPu , GPm , GPd ;
	double r , mu , prior , maxprior , sigmabar, minprior , dt=0.0 , S , t;
	double Ti , Ki , Ci ;
  	double Ggama, gama; // derive seconde de V %S  (Vss)
	double Ttemp , temp = 0.0 , Gtemp = 0.0 , temp0 = 0.0 ;
	char *flag_payoff ;


	flag_payoff = (char *)malloc(sizeof(char));
	tdim = arbre.tdim ;

	n = tdim ;// le nombre de date dans V, (n+1)^2=taille de V

	costFunction(lam);
	r = param.r ;
	mu = param.mu ;
	prior = param.prior ;
	sigmabar = param.sigmabar ;
	minprior = param.prior_min ;
	maxprior = param.prior_max ;
			
	Ttemp = SOption.maturite[0] ;
	
	for( k = 0 ; k <SOption.size ; k++)
	{
		Ki = SOption.strik[k] ;
		Ci= SOption.prix[k] ;
		Ti = SOption.maturite[k] ;
		if (Ti==Ttemp) rank = (int)indice.val[m];
		else Ttemp = Ti, m++,rank = (int) indice.val[m];
            
		*flag_payoff = SOption.payoff[k] ;
		t = Ti ;

		for( i =0 ; i < (tdim+1)*(tdim+1) ; i++)
    	{
    		Gsource[i] = 0.0 ;
    	}

		define_G_source( Gsource, Ki , flag_payoff , arbre , rank ) ;
		for ( i = 0 ; i< 2 ; i++ )
		{
			for( j = 0 ; j < (tdim+1)*(tdim+1) ; j++) GV.val[i][j] = 0.0 ;
		}
				
		//valeur terminale de GV a la date T, la plus grande des maturites
		i = rank ;
		for ( j = i*i ; j < (i+1)*(i+1) ; j++)
		{
			GV.val[0][j] = Gsource[j] ;
		}
				
		for ( i = rank-1 ; i >= 0 ; i--)
    	{
    		dt = arbre.pasTemps[i*i] ;
    		t-= dt ;
       		for ( j = i*i ; j < (i+1)*(i+1) ; j++)
       		{
         		// calcul du parametre p(i,j) , Gp(i,j) a partir de phi , phi' et de V a la date i+1 et au niveau (j,j+1,j+2)
        	  	//p(i,j)=(phi(exp(-0.5rt S^2 Vss)/(sigmabar^2 exp(-rt) 0.5 S^2 Vss)
      	    	S = arbre.val[j] ;
         		temp0 = sigmabar*sqrt(dt)/2.0 ;
         		gama = (( 1 - temp0 )*V.val[0][2*i+j+3]+( 1.0 + temp0 )*V.val[0][2*i+j+1]- 2.0*V.val[0][2*i+j+2] )/(4*temp0*temp0);
         		Ggama =(( 1 - temp0 )*GV.val[0][2*i+j+3]+( 1.0 + temp0 )*GV.val[0][2*i+j+1]- 2.0*GV.val[0][2*i+j+2] )/(4*temp0*temp0);
	          	temp = exp(-r*t)*Ggama ;
				
    	  	    //Gp = GPHI( temp , prior , minprior , maxprior )/ ( sigmabar*sigmabar ) ;
         		Gp = GPHI( exp(-r*t)*gama , prior , minprior , maxprior )/ ( sigmabar*sigmabar ) ;
         		//printf("gama %f\n", gama);
		        GV.val[1][j] = Gp ;
                    
        	  	//calcul des proba Pu(i,j) , Pm(i,j), Pd(i,j) et de  GPu(i,j) , GPm(i,j), GPd(i,j)
		        //Pu=0.5p(1-0.5sigmabar sqrt(dt) ) +0.5 mu dt/sigmabar
     		    //Pm=1-p
         		//Pd= 0.5  p(1+0.5sigmabar sqrt(dt) ) -0.5mu sqrt(dt)/sigmabar
            	// meme chose pour GP avec Gp
	            GPu = 0.5*Gp*(1-0.5*sigmabar* sqrt(dt)) +0.5*mu*sqrt(dt)/sigmabar ;
     		    GPm=1-Gp;
         		GPd= 0.5*Gp*(1+0.5*sigmabar*sqrt(dt)) -0.5*mu*sqrt(dt)/sigmabar ;
				if(GPu<=0 || GPu>=1  || GPm>=1 || GPm<=0 || GPd>=1|| GPd<=0  ) printf ("\n\n\nPROBA NEGATIVE dans GV GPU %f GPm %f GPd %f \n\n",GPu,GPm,GPd);
    	
    	        //calcul de GV(i,j)=GPu(i,j)GV(i+1,j+2)+GPm(i,j)GV(i+1,j+1)+GPd(i,j)GV(i+1,j)+Gsource(i,j)
 		        GV.val[0][j] = GPu* GV.val[0][2*i+j+3]+GPm*GV.val[0][2*i+j+2]+GPd*GV.val[0][2*i+j+1] ;
	    		GV.val[0][j] = GV.val[0][j]*exp(-r*dt)+ Gsource[j];

      		}
    	}
					
		grad[k] = GV.val[0][0]-Ci;
	}

	free(flag_payoff);

}

double price( double K , double S, double T , Parametre param)
{

   	Trinomial_tree Tree ;
   	Tri_tree_prob_space Vprix ;
   	Vecteur diffMat;
   	Vecteur indis;
   	int i,j=0 ,n, dim = 1;
   	double s , res , Pm,Pu,Pd,temp0,t,p,dt;
   	char *F;

   	indis.size = dim ;
   	indis.val = (double *) malloc(dim*sizeof(double));
   	indis.val[0] = 0 ;
   	diffMat.size = dim ;
   	diffMat.val = (double *) malloc(dim*sizeof(double));
   	diffMat.val[0] = T ;
   	dt=T/param.N;
   	param.dt=dt;
   	Tree = getPriorTree(diffMat,&indis , param);
   	//Tree = getPriorTree(param);

   	Vprix = define_3_tree_prob_space( Tree.tdim , 0.0 );

   	n=i=Tree.tdim;
   	F=(char*)malloc(sizeof(char));
   	*F= SOption.payoff[0] ;
   	printf("dans price\n");

   	printf("S %f K %f T %f payoff %c\n",S,K,T,SOption.payoff[0]);
   	t=T;
   	for ( j = n*n ; j < (n+1)*(n+1) ; j++)
  	{
	    s = Tree.val[j];
        Vprix.val[0][j] = Gi(s,K,F);
	}

  	for ( i = n-1 ; i >= 0 ; i--)
    {
   	 	dt = Tree.pasTemps[i*i] ;
		t-= dt ;

      	for ( j = i*i ; j < (i+1)*(i+1) ; j++)
        {
   			s = Tree.val[j] ;
   			temp0 = param.sigmabar*sqrt(dt)/2.0 ;

   			p = param.prior*param.prior/(param.sigmabar*param.sigmabar);

   			//p = 1 ;
		    Pu = 0.5*p*(1-temp0) +0.5*param.mu*sqrt(dt)/param.sigmabar ;
            Pm = 1-p ;
            Pd = 0.5*p*(1+ temp0) -0.5*param.mu*sqrt(dt)/param.sigmabar ;
	    	Vprix.val[0][j] = Pu* Vprix.val[0][2*i+j+3]+Pm*Vprix.val[0][2*i+j+2]+Pd*Vprix.val[0][2*i+j+1] ;
   			Vprix.val[0][j] = Vprix.val[0][j];

      	}
    }

   	res = Vprix.val[0][0]*exp(-param.r*T) ;
   	free_Tri_tree_prob_space(Vprix) ;
    free_Trinomial_tree(Tree);
	free(F) ;
    free_Vecteur(indis);
    free_Vecteur(diffMat);
	return(res) ;
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
