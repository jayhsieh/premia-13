# include "maths.h" 

 unsigned long   u=12345;

unsigned long getrandmax(){
return 2147483647;
}


double rnd(){
u=(16807*u)%getrandmax() ;
return u;
}

/**simulate a uniform random between 0 and 1**/

double unif(){
double a;
a= (rnd()/getrandmax());
return a;
}
/** simulate a gaussian random**/

double gaussian(){
double b;
double pi=3.14159265;
b=sqrt(-2*log(unif())) *cos(2*pi*unif());
return b;
}

