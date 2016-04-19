
#ifndef common_h
#define common_h

#include <vector>

long double cumnorm(long double&);
long double deval (long double&, long double *, int&);
long double neval (long double&, long double *, int&);
long double p1evll(long double&, void *PP, int&);
long double polevll(long double&, void *PP, int&);
long double erfcel(long double);
long double dnorm(long double);

//fonction Hh-1
long double Hh(const long double&);

//fonction Hh0
long double Hh0(const long double&);

//fonctions qui retourne le vecteur des n premiers fonctions Hhn
std::vector<long double> Hhn(const long double&,const int&);

//fonction I-1
long double I(const long double&,const long double&,const long double&,const long double&);

//fonction I0
long double I0(const long double&,const long double&,const long double&,const long double&);

//fonction qui retourne le vecteur des n premiers fonctions In
std::vector<long double> In(const long double&,const long double&,const long double&,const long double&,const int&);

//fonction qui retourne le vecteur des n premiers fonctions derivés de In % a sa première variable
std::vector<long double> dIn(const long double& ,const long double& ,const long double& ,const long double& ,const int&);

//fonction factorielle
unsigned long long int fact_dia(const int &);

//coefficients binomiaux
unsigned long long int bin_dia(const int&,const int&);


#endif
