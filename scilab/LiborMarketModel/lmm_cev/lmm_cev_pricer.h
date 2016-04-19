#ifndef _LMM_CEV_
#define _LMM_CEV_

#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <time.h>


/*Pricing function*/
//int_pts= MC or number of discretiz fore auxiliary functions
void cev_price(size_t type_model,size_t type_product,size_t type_pricing, 
		 size_t type_scheme, double alpha, double fzero, double H,
	       double expiry,double swp_lenght,int int_pts,double* out);

/*auxiliary functions to make pricer self-consistent, stand-alone*/
// Random generation
int function_n(double t);
double u_random();
double simul_normale();
//miscellanea
double function_fi(double x,int i,double alpha);
double function_delta( int k);
double function_lambda_j(double t,int k);
double function_lambda_k_carre(double t);
double price_zero_coupon_bond(double t,int k,double fzero);
double function_B_s(double t,int debut,int fin,double fzero);
double function_R(double t,int debut,int fin,double fzero);
double function_dr_df(double t,int j,int debut,int fin,double fzero);
double function_wj(double t, int j,double alpha,int debut,int fin,int i,double fzero);
double function_somme(double t,int debut,int fin,int i,double alpha,double fzero);
double somme_integral(double a,double b,int debut,int fin,int i,double alpha,double fzero,int n);
double calcul_integral(double (*f) (double),double a, double b, int n);
int function_a(int k, double t,double alpha,double H,int n);
double function_b(double alpha);
double function_d(double principal_teta,double t,double vs,double alpha);
double function_petitf(double t,double vs,double alpha,int debut,int fin,double fzero);
double function_gplus(double t,double principal_teta,double vs,int debut,int fin,double fzero);
double function_gmoins(double t,double principal_teta,double vs,int debut,int fin,double fzero);
double function_c(int k,double t,double alpha,double fzero,int n);
double function_xplus(double t,int k,double alpha,double H,double fzero,int n);
double function_xmoins(double t, int k,double alpha,double H,double fzero,int n);
//ChiSquared simulation
double distrib_normale_positive(double x);
double distrib_normale_negative(double x);
double distrib_normale(double x);
double alngam(double xvalue)  ;
double distrib_chi2(double x,double theta/* degrés */,double f/* decentrage */);
//Low level pricing functs
double prix_caplet(double t,int k,double alpha,double H,double fzero,int n);
double prix_swaption(double t,double principal_teta,double vs,double alpha,int debut,int fin,double fzero);
double Black(int k,double sigma,double H,double fzero);
double Black_prime(int k,double sigma,double H,double fzero);
double implied_vol(double sigma_n,int k,double prix_caplet,double H,double fzero);
void function_F(int e,double t[],double final[],double fzero,double delta_i,double alpha,int type_model,int type_scheme);
//High level pricing function
void function_F(int e,double t[],double final[],double fzero,double delta_i,double alpha,int type_model,int type_scheme);
void function_F_swap(int debut,int fin,double t[],double final[],double fzero,double delta_i,double alpha,int type_model,int type_scheme,double finalcor[]);
double function_B_en_maturite(double o,double F[]);
#endif //_LMM_CEV_


