#ifndef TREE
#define TREE

#define _b1 0.2316419
#define _b2 0.31938153
#define _b3 -0.356563782
#define _b4 1.781477937
#define _b5 -1.821255978
#define _b6 1.330274429
#define _pi  3.1415926535
#define _sqrtdeuxpi 2.5066282746
#define _PrecisionVol 0.0001

typedef struct {
  				int tdim;//si n=2 le nombre de noeuds=9
				double *val ;
				double *pasTemps ;
				double *proba;
				}Trinomial_tree;
											
typedef struct {
    			int tdim;//si n=2 le nombre de noeuds=9
				double **val ;
              	}Tri_tree_prob_space;


typedef struct {
                double *maturite ;
                double *strik ;
                double *prix ;
                char *payoff ;
                int size ;
                } Vect_option ;

typedef struct {
                double r ;
                double mu ;
                double dividende ;
                double T ;
                double S0 ;
                double N ;
                double prior ;
                double dt ;
                double sigmabar ;
                double prior_min ;
                double prior_max ;
                } Parametre;
typedef struct {
                int size ;
                double *val ;
                 } Vecteur ;


Trinomial_tree define_3_tree(/*! nb d'interval*/ int tdim , double defaultVal );

Tri_tree_prob_space define_3_tree_prob_space(int tdim , double defaultVal );
Vecteur define_vecteur( int ) ;

Vect_option define_option( const int dim ) ;
Vect_option load_options( char *file_title);

int comp_nums(const double *num1, const double *num2);

Vecteur getSortedMaturity( Vect_option );
Vecteur getDiffMaturity( Vect_option );

Vect_option Sort( Vect_option option );
/**************************************************************************
*  Gi retourne (S-K)+ pour un  Call ou un Put, F=C ou F=P *
****************************************************************************/
double Gi(double S,double K,char *F);

void affiche_tree( Trinomial_tree ) ;
void afficheVectTree( double * vecteur , int lengthTree ) ;

void free_Trinomial_tree(Trinomial_tree);
void free_Tri_tree_prob_space(Tri_tree_prob_space);
void free_Vect_option(Vect_option);
void free_Vecteur(Vecteur);


//fonction de répartition de la loi normale N(0,1)
double Rep_Normale(double d) ;
//double distrib(double d);
double Max(double a,double b) ;
double Min(double a,double b) ;
double BS_d1(double S,double VolImpl,double r,double K,double T,double q) ;

double BS_d2(double S,double VolImpl,double r,double K,double T,double q) ;

double BS(char *C,double S,double VolImpl,double r,double K,double T,double q);

double ImpliedVol(double S_0, double r, double q, char *optionType, double t_0, double K, double T, double V, double tol);
double ImpliedVolNewton(double S_0, double r, double q, char *optionType, double t_0, double K, double T, double V, double tol);
int pgcd(int a , int c);
double getdt(Vecteur differenteMaturite,int precison);

#endif
