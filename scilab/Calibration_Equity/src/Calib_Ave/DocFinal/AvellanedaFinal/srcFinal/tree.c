#include "tree.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double Gi(double S,double K,char *F)
{
    double res=0.0;
    //double S0=param.S0;
    if (*F=='C')
    {
				res = (S-K)*(S>=K);
    }
    if (*F=='P')
    {
				res = (K-S)*(S<=K);
    }
    return res;
}

Trinomial_tree define_3_tree(/*! nb d'interval*/ int tdim , double defaultVal )
{
	int i ;
	Trinomial_tree tree_nom;
	tree_nom.tdim = tdim;
	tree_nom.val = (double *)calloc((tree_nom.tdim+1)*(tree_nom.tdim+1),sizeof(double));
	tree_nom.pasTemps = (double *)calloc((tree_nom.tdim+1)*(tree_nom.tdim+1),sizeof(double));
	tree_nom.proba=(double *)calloc((tree_nom.tdim+1)*(tree_nom.tdim+1),sizeof(double));
    for( i =0 ; i < (tree_nom.tdim+1)*(tree_nom.tdim+1) ; i++)
    {
    	tree_nom.val[i] = defaultVal ;
    	tree_nom.pasTemps[i] = 0.0 ;
    	tree_nom.proba[i] = 0.0 ;
    }

	return (tree_nom);
}			

Tri_tree_prob_space define_3_tree_prob_space( int tdim , double defaultVal )
{
		int i , j ;
		Tri_tree_prob_space tree_nom;
		
		tree_nom.tdim=tdim;
		tree_nom.val=(double **) malloc((2)*sizeof(double *));

		for ( i = 0 ; i<= 1 ; i++ )
		{
			tree_nom.val[i] = (double *) malloc(((tree_nom.tdim+1)*(tree_nom.tdim+1))*sizeof(double));
			for( j = 0 ; j < (tdim+1)*(tdim+1) ; j++) tree_nom.val[i][j] = defaultVal ;
		}

		return (tree_nom);

}

Vect_option define_option( const int dim )
{
   	Vect_option temp_option ;
   	temp_option.size = dim ;//si on a n option dim =n
   	temp_option.maturite =(double *) malloc(dim*sizeof(double)) ;
   	temp_option.prix =(double *) malloc(dim*sizeof(double)) ;
   	temp_option.strik =(double *) malloc(dim*sizeof(double)) ;
   	temp_option.payoff =(char *) malloc(dim*sizeof(char)) ;
   	printf("option size dans define %d dim %d\n",temp_option.size,dim);
   	return (temp_option) ;
    free_Vect_option ( temp_option ) ;
}

Vecteur define_vecteur( int dim )
{
    Vecteur vect ;
    vect.size = dim ;
    vect.val = (double *) malloc(dim*sizeof(double));

    return vect ;
}

Vect_option load_options(char *file_title )
{
    double Ki,prix,Ti ;
	//char payoff ;
	char flag;
	int i = 0 , j = 0 ;
	char ch_line[200];
	Vect_option option;
	FILE *fi ;

	while((fi=fopen(file_title,"r"))==NULL)
		printf("impossible d'ouvrir le fichier \n");

    // get number of line
   	while( fgets(ch_line,200,fi)!= NULL ) j++;
   		fclose(fi);

    printf("taille du fichier %d\n",j);
    //fin get number of line

   	while((fi=fopen(file_title,"r"))==NULL)
		printf("impossible d'ouvrir le fichier \n");

   	option.size = j ;//si on a n option dim =n
   	option.maturite =(double *) malloc(j*sizeof(double)) ;
	option.prix =(double *) malloc(j*sizeof(double)) ;
	option.strik =(double *) malloc(j*sizeof(double)) ;
	option.payoff =(char *) malloc(j*sizeof(char)) ;

   	//fgets (ch_line,80, *file_title);
   	//dans le cas ou on a une indication sur le payoff(Call, Put)
   	//char flag, opt.pay[i]=flag
	while(fscanf(fi,"%lf %lf %lf %c\n", &Ki,&Ti,&prix,&flag)!=EOF)
	//while(fscanf(fi,"%lf %lf %lf", &Ki,&Ti,&prix)!=EOF)
	{
		//printf("%f %f %f \n",Ki,Ti,prix);
		option.maturite[i] = Ti ;
		option.prix[i] = prix ;
		option.strik[i] = Ki ;
		//flag ='C' ;
		option.payoff[i]=flag ;
		i++ ;
	}

	//free(flag);
	fclose(fi);
	printf("succees: fin load, taille option %d\n ",option.size);
	return option ;

}

/************************************************************************
 *
 * comp_nums: Compare two numbers.
 *
 ************************************************************************/

int comp_nums(const double *num1, const double *num2)
{
if ( *num1 <  *num2 )
  {
    return -1;
  }
  else if ( *num1 == *num2 )
  {
    return  0;
  }
  else if (*num1 >  *num2) return  1;
  {
    return  1;
  }
	/*if (*num1 <  *num2) return -1;
	if (*num1 == *num2) return  0;
	if (*num1 >  *num2) return  1;*/
}

/************************************************************************
 *
 * getSortedMaturity: renvoie un pointeur sur le vecteur des maturites
 * dans l'ordre croissant.
 *
 ************************************************************************/

Vecteur getSortedMaturity( Vect_option option )
{
	Vecteur res ;
    int how_many = option.size ;
    int i  ;
    res = define_vecteur( how_many ) ;

    for( i = 0; i<how_many;i++)
    {
    	res.val[i]=option.maturite[i] ;
    }

    qsort(
        res.val, 			/* Pointer to elements		*/
        how_many, 			/* Number of elements		*/
        sizeof(double),  			/* size of one element.		*/
        (void *)comp_nums		/* Pointer to comparison function */);
	return (res) ;

}

/************************************************************************
 *
 * getDiffMaturity: renvoie un pointeur sur le vecteur DIFFERENTES
 *  maturites dans l'ordre croissant.
 *
 ************************************************************************/
Vecteur getDiffMaturity( Vect_option option )
{
	Vecteur res , diffMatu ;
    int how_many = option.size ;
    int i , j ;

    res = define_vecteur( how_many ) ;

    for( i = 0; i<how_many;i++)
    {
        res.val[i]=option.maturite[i] ;
    }

    qsort(res.val, how_many, sizeof(double), (void *)comp_nums);
    j = 1 ;
    for(  i = 0; i<how_many-1;i++)
    {
	    if (res.val[i+1]!=res.val[i]) j++ ;
    }

    diffMatu = define_vecteur(j) ;

    j = 0 ;
	diffMatu.val[j] = res.val[j] ;
    
    for(  i = 0; i<how_many;i++)
    {
    	if (res.val[i]!=diffMatu.val[j])
        {
        	diffMatu.val[j+1] = res.val[i] ;
           	j++ ;
        }
	}
    free_Vecteur(res) ;
    return ( diffMatu) ;
	
}

/************************************************************************
 *
 * Sort: renvoie une structure vect_option copie d'une option ordonnee selon les maturites croissantes
 *  la structure passee en parametre n'est pas modifiee.
 *
 ************************************************************************/
Vect_option Sort( Vect_option option )
{
	int i, j, inc , k;
    Vect_option tempoption ;
    int n = option.size;
    double v,w,x;
    char z;

    tempoption=define_option(n);
    for( k = 0 ; k <tempoption.size ; k++)
    {
        tempoption.maturite[k]=option.maturite[k];
        tempoption.payoff[k] = option.payoff[k];
        tempoption.prix[k] = option.prix[k];
        tempoption.strik[k] = option.strik[k];
    }

    inc = 1;

    do
    {
        inc *= 3;
        inc++;
    }
    while (inc <= n);

    do
    {
        inc /= 3;

        for (i=inc+1; i<=n; i++)
        {
            //v = array.Elt(i-1);
            v = tempoption.maturite[i-1] ;
            w = tempoption.strik[i-1] ;
            x = tempoption.prix[i-1] ;
            z = tempoption.payoff[i-1] ;

            j=i;

            while (tempoption.maturite[j-inc-1] > v)
            {
                tempoption.maturite[j-1] = tempoption.maturite[j-inc-1];
                tempoption.strik[j-1] = tempoption.strik[j-inc-1];
                tempoption.prix[j-1] = tempoption.prix[j-inc-1];
                tempoption.payoff[j-1] = tempoption.payoff[j-inc-1];

                j -= inc;

                if (j <= inc) break;
            }
             tempoption.maturite[j-1] = v ;
             tempoption.strik[j-1] = w ;
             tempoption.prix[j-1] = x ;
             tempoption.payoff[j-1] = z ;
        }
    }
    while (inc>1);
    printf("succeed end of sort . Size of sorted option %d\n",tempoption.size);
    return(tempoption);
}

void free_Trinomial_tree(Trinomial_tree tree)
{
    free(tree.val);
    free(tree.pasTemps);
    free(tree.proba);
}

void free_Tri_tree_prob_space(Tri_tree_prob_space tree)
{
    free(tree.val[0]);
    free(tree.val[1]);
    free(tree.val);
}

void free_Vect_option(Vect_option option)
{
    free(option.maturite);
    free(option.payoff);
    free(option.prix);
    free(option.strik);
}

void free_Vecteur(Vecteur vect)
{
    free(vect.val);
}


double Rep_Normale(double d)
{
    double x, y, z ;
    if (d >= 0)
	{
		x = 1 / (1 + _b1 * d) ;
    y = exp(-d * d / 2.) / _sqrtdeuxpi ;
    z = x * (_b2 + x * (_b3 + x * (_b4 + x * (_b5 + x * _b6)))) ;
    z = 1 - y * z ;
	}
    else
    {
		x = 1 / (1 - _b1 * d) ;
    y = exp(-d * d / 2) / _sqrtdeuxpi ;
		z = x * (_b2 + x * (_b3 + x * (_b4 + x * (_b5 + x * _b6)))) ;
		z = y * z ;
    }
	return (z) ;
}


double Max(double a,double b){ return (b>a)?b:a;}

double Min(double a,double b){ return (b>a)?a:b;}

double BS_d1(double S, double VolImpl, double r, double K, double T , double q)
{
	return ((log(S/K)+r*T-q*T+VolImpl*VolImpl*T/2.) / VolImpl / sqrt(T)) ;
}

double BS_d2(double S, double VolImpl, double r, double K, double T, double q)
{
	double resu ;
	resu =  BS_d1(S,VolImpl,r,K,T,q) - VolImpl * sqrt(T) ;
	return resu ;
}

// formule de BS tout en vrai valeur (pas en %)
// int C : C=1 call sinon put
double BS(char *C,double S,double VolImpl,double r,double K,double T , double q)
{
	double resu ,d1 ,d2 ,Nd1 , Nd2;
	
	if ( *C =='C')
	{
		if (T==0. || VolImpl == 0.)
			resu = Max(S-K,0.) ;
		else
		{
			
			d1 = BS_d1(S,VolImpl,r,K,T,q) ;
			d2 = d1 - VolImpl * sqrt(T) ;
			
			Nd1 = Rep_Normale(d1);
			Nd2 = Rep_Normale(d2);
			
			resu = S*exp(-q*T)*Nd1-K*Nd2*exp(-r*T) ;
		}
	}
	else
	{
		
		if (T==0.||VolImpl==0.)
			resu = Max(K-S,0.) ;
		else
		{
			d1 = BS_d1(S,VolImpl,r,K,T,q) ;
			d2 = d1 - VolImpl * sqrt(T) ;
			Nd1 = Rep_Normale(-d1);
			Nd2 = Rep_Normale(-d2);
			resu = K*exp(-r*T)*Nd2-S*exp(-q*T)*Nd1 ;
		}
	}

	return resu ;
}


void affiche_tree( Trinomial_tree Tree )
{
		int i , j ;
		for ( i = 0 ; i < Tree.tdim + 1 ; i++)
    {
      for ( j = i*i ; j < (i+1)*(i+1) ; j++)
      {
        printf(" %f",Tree.val[j]) ;
      }
      printf("\n");
    }
}

void afficheVectTree( double * vecteur , int lengthTree )
{
		int i , j ;
		for ( i = 0 ; i <  lengthTree + 1 ; i++)
    {
      for ( j = i*i ; j < (i+1)*(i+1) ; j++)
      {
        printf("%0.3f ",vecteur[j]) ;
      }
      printf("\n");
    }
}



double ImpliedVol(double S_0, double r, double q, char *optionType, double t_0, double K, double T, double V, double tol)
{
    double sigma_min,sigma_max;
    double sigma;
    double price;
    int IterMax=100,iter=0;

    sigma_min = 0;
    sigma_max = 10;
    sigma = 0.5*(sigma_max+sigma_min);

    if ( *optionType=='C' && V < S_0*exp(-q*T)-K*exp(-r*T)) return 0.0;/*impossible pour un prix BS d'un Call*/
    if ( *optionType=='P' && V < -S_0*exp(-q*T)+K*exp(-r*T)) return 0.0;/*impossible pour un prix BS d'un Putt*/
	price = BS(optionType,S_0,sigma , r , K , T , q) ;

	while (fabs(price-V)>tol && iter<IterMax)
  	{
    	if (price-V > 0) sigma_max = sigma;
	    else	sigma_min = sigma;
  		sigma = 0.5*(sigma_max+sigma_min);
	   	price = BS(optionType,S_0,sigma , r , K , T , q);
   		iter++;
  	}

	if (iter<IterMax ) {return sigma;}
	else
	{
		printf ("problem dans vol implicite pour S0 %f, r %f, q %f, t0 %f, K %f, T %f, Prix %f, Tol %f, ImpVol %f",S_0,r,q,t_0,K,T,V,tol,sigma);
		return -1e40;
	}

}

double ImpliedVolNewton(double S_0, double r, double q, char *optionType, double t_0, double K, double T, double V, double tol)
{
    double sigma_min = 1e-5,sigma_max;
    double sigma;
    double price= BS(optionType,S_0,sigma_min , r , K , T , q);
    double t_sqrt = sqrt(T);
    double diff,d1,d2,vega;
    int i,IterMax=100,iter=0;

	if (V < S_0*exp(-q*T)-K*exp(-r*T)) return 0.0;/*impossible pour un prix BS*/

    sigma_max = 10;
    sigma =(V/S_0)/(0.398*t_sqrt);
	sigma =5;
	for (i=0;i<IterMax;i++)
	{
    	price = BS(optionType,S_0,sigma,r,K,T,q);
	    diff = V - price;

    	if (fabs(diff)<tol) return sigma;

	    d1 = (log(S_0/K)+r*T)/(sigma*t_sqrt) + 0.5*sigma*t_sqrt;
		d1= BS_d1(S_0, sigma, r, K, T, q);
    	d2= BS_d2(S_0, sigma, r, K, T, q);
	    vega = -S_0 *exp(-q*T)*d2*exp(-d1*d1/2)/sigma/_sqrtdeuxpi +K*exp(-r*T)*d1*exp(-d2*d2/2)/sigma/_sqrtdeuxpi;
    	sigma = sigma + diff/vega;
  	};
  	return -99e10;  // something screwy happened, should throw exception

}


int pgcd(int a , int b)
{
	int res=0;
	if(a==b) res=a;
	else
	{
		if ( a > b ) res = pgcd(a-b,a);
		else res = pgcd(a,b-a);
	}
	return res;
}

double getdt(Vecteur differenteMaturite, int precision)
{
	Vecteur vect;
	int i,m=0,max=0,res=0;


	vect.size = differenteMaturite.size ;
	vect.val = (double *) malloc(vect.size*sizeof(double));
	for (i=0;i<vect.size;i++)
	{
		vect.val[i] = ceil(pow(10,precision)*differenteMaturite.val[i]);
	}

	res=(int)vect.val[vect.size-1];
	for(i=vect.size-1;i>0;i--) res=pgcd((int)vect.val[i-1],res);
	free(vect.val);
 	return res/(double)pow(10,precision);
}
