#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_cdf.h"
#include "pnl/pnl_integration.h"
#include "pnl/pnl_root.h"
#include "swaption_hw2d.h"
#include "swaption_black.h"

static double xt,yt,a,b,sigma,eta,rho;
static double tau,Nominal,K;
static double omega;
static double option_mat,swap_mat;
static double mu_x,mu_y,sigma_x,sigma_y,rho_xy;
static double critical_y;
static int  nb_payement;

static PnlVect* Ci;
static PnlVect* log_Ci;
static PnlVect* log_A_func_T_ti;
static PnlVect* B_func_a_T_ti;
static PnlVect* B_func_b_T_ti;

static double V_func(double t,double T)
{
    return SQR(sigma/a)*(T-t+2./a*exp(-a*(T-t))-1./(2.*a)*exp(-2*a*(T-t))-3./(2*a))+
           SQR(eta/b)*(T-t+2./b*exp(-b*(T-t))-1./(2.*b)*exp(-2*b*(T-t))-3./(2*b))+
           2.*rho*sigma*eta/(a*b)*(T-t+(exp(-a*(T-t))-1.)/a+(exp(-b*(T-t))-1.)/b-(exp(-(a+b)*(T-t))-1.)/(a+b));
}

static double log_A_func(ZCMarketData* ZCMarket, double t, double T)
{
    double VtT,V0T,V0t;
    double P0_t,P0_T;

    VtT= V_func(t,T);
    V0T= V_func(0,T);
    V0t= V_func(0,t);

    P0_T = BondPrice(T, ZCMarket);
    P0_t = BondPrice(t, ZCMarket);

    return log(P0_T/P0_t) + 0.5*(VtT-V0T+V0t);
}

static double B_func(double z,double t,double T)
{
    if (z==0) return (T-t);
    return (1./z)*(1.-exp(-z*(T-t)));
}

static double MxT(double s,double t,double T)
{
    return (SQR(sigma/a)+rho*sigma*eta/(a*b))*(1.-exp(-a*(t-s)))
           -0.5*SQR(sigma/a)*(exp(-a*(T-t))-exp(-a*(T+t-2*s)))
           -(rho*eta*sigma/(b*(a+b)))*(exp(-b*(T-t))-exp(-b*T-a*t+(a+b)*s));
}

static double MyT(double s, double t, double T)
{
    return  (SQR(eta/b)+rho*sigma*eta/(a*b))*(1.-exp(-b*(t-s)))
            -0.5*SQR(eta/b)*(exp(-b*(T-t))-exp(-b*(T+t-2*s)))
            -(rho*eta*sigma/(a*(a+b)))*(exp(-a*(T-t))-exp(-a*T-b*t+(a+b)*s));
}

static double h1(double x)
{
    return (critical_y-mu_y)/(sigma_y*sqrt(1.-SQR(rho_xy)))-rho_xy*(x-mu_x)/(sigma_x*sqrt(1.-SQR(rho_xy)));
}

static double log_lamda_func(ZCMarketData* ZCMarket, double x,double ti,int i)
{
    return  GET(log_Ci,i) + GET(log_A_func_T_ti,i) - GET(B_func_a_T_ti,i)*x;
}

static double ki_func(double x, int i)
{
    return -GET(B_func_b_T_ti,i)*(mu_y-0.5*(1.-SQR(rho_xy))*SQR(sigma_y)*GET(B_func_b_T_ti,i)+rho_xy*sigma_y*(x-mu_x)/sigma_x);
}

/*Computation of Critical Y*/
static double phiY(double y, void *p)
{
    int i;
    double x=*((double*)p), sum, tmp, phi_y;
    sum=0.;
    for (i=0; i<nb_payement; i++)
    {
        tmp = GET(Ci, i)*exp(GET(log_A_func_T_ti,i) - GET(B_func_a_T_ti,i)*x - GET(B_func_b_T_ti,i)*y);
        if (tmp>2.) return 10.;
        sum += tmp;
    }
    phi_y = sum-1.;

    return phi_y;
}


static void bracket_phi(double *y_1, double *y_2, double x)
{
    double y, y_step=10, phi_value;

    *y_1 = (GET(log_A_func_T_ti, nb_payement-1) - GET(B_func_a_T_ti, nb_payement-1)*x)/GET(B_func_b_T_ti, nb_payement-1);
    if (*y_1>0) *y_1 *= 0.5;
    else *y_1 *= 2;
    phi_value = phiY(*y_1, &x);

    y=*y_1;
    y_step = MAX(100., fabs(y));
    while (phi_value<0)
    {
        y = y-y_step;
        phi_value = phiY(y, &x);
        *y_1 = y;
    }

    y=*y_1;
    while (phi_value>0)
    {
        y = y+y_step;
        phi_value = phiY(y, &x);

        *y_2 = y;
    }

}

static double Critical_Y_Bisection(ZCMarketData* ZCMarket, double x)
{
    int N_max=1000;
    double epsabs=1e-10, epsrel=1e-10;
    double y_0=0., y_1=0., y_2=0.;
    PnlFunc func;

    func.function = phiY;
    func.params = &x;

    bracket_phi(&y_1, &y_2, x);
    pnl_root_bisection(&func, y_1, y_2, epsrel, epsabs, N_max, &y_0);

    return y_0;
}
// Function to integrate from ]-inf, inf[
static double integrand_fun(double x, void* ZCMarket)
{
    double d1,d2,sum, tmp1, tmp2, res;
    double ti;
    int i;

    critical_y=Critical_Y_Bisection((ZCMarketData*)ZCMarket,x);
    d1=-omega*h1(x);
    sum=0.;
    ti = option_mat;
    tmp1 = -0.5*SQR((x-mu_x)/sigma_x);
    tmp2 = -omega*sigma_y*sqrt(1.-SQR(rho_xy));
    for (i=0; i<nb_payement; i++)
    {
        ti += tau;
        d2 = d1 + GET(B_func_b_T_ti, i)*tmp2;
        sum += exp(tmp1+ki_func(x,i) + log_lamda_func((ZCMarketData*)ZCMarket,x,ti,i))*cdf_nor(d2);
    }

    res = exp(tmp1)*cdf_nor(d1) - sum;

    return res;
}

/*Payer Swaption momega=1, Receiver momega=-1*/
static double SWAPTION_g2(ZCMarketData* ZCMarket, double momega, double moption_mat, double mtau, double mswap_mat, double mNominal, double mK, double mxt, double myt, double ma, double mb, double msigma, double meta, double mrho)
{
    double sum;
    double Tim, P0_Ti, P0_opmat, P0_swapmat, S0_Ti;
    int i, neval;
    double integral, abserr=0.;
    double epsabs=1e-7, epsrel=1e-7;
    int limit=0;
    PnlFunc func;

    omega=momega;
    option_mat=moption_mat;
    tau=mtau;
    swap_mat=mswap_mat;
    Nominal=mNominal;
    K=mK;
    xt=mxt;
    yt=myt;
    a=ma;
    b=mb;
    sigma=msigma;
    eta=meta;
    rho=mrho;

    nb_payement=(int)((swap_mat-option_mat)/tau);
    Ci =  pnl_vect_create_from_double(nb_payement, K*tau);
    LET(Ci, nb_payement-1) = 1+K*tau;

    P0_opmat = BondPrice(option_mat, ZCMarket);
    P0_swapmat = BondPrice(swap_mat, ZCMarket);

    log_Ci =  pnl_vect_create(nb_payement);
    log_A_func_T_ti =  pnl_vect_create(nb_payement);
    B_func_a_T_ti =  pnl_vect_create(nb_payement);
    B_func_b_T_ti =  pnl_vect_create(nb_payement);

    sum=0.;
    Tim = option_mat;
    for (i=0;i<nb_payement;i++)
    {
        Tim += tau;
        P0_Ti = BondPrice(Tim, ZCMarket);
        sum += tau*P0_Ti;

        LET(log_Ci, i) = log(GET(Ci, i));
        LET(log_A_func_T_ti, i) = log_A_func(ZCMarket, option_mat, Tim);
        LET(B_func_a_T_ti, i) = B_func(a, option_mat, Tim);
        LET(B_func_b_T_ti, i) = B_func(b, option_mat, Tim);
    }

    /*Compute Swap Rate*/
    S0_Ti = (P0_opmat-P0_swapmat)/sum;

    /*Integral computation*/
    mu_x=-MxT(0,option_mat,option_mat);
    mu_y=-MyT(0,option_mat,option_mat);

    sigma_x=sigma*sqrt((1.-exp(-2.*a*(option_mat)))/(2.*a));
    sigma_y=eta*sqrt((1.-exp(-2.*b*(option_mat)))/(2.*b));
    rho_xy=rho*eta*sigma/((a+b)*sigma_x*sigma_y)*(1.-exp(-(a+b)*option_mat));

    func.function = integrand_fun;
    func.params = ZCMarket;

    pnl_integration_qag(&func, PNL_NEGINF, PNL_POSINF, epsabs, epsrel, limit, &integral, &abserr, &neval);
    integral = integral/(sigma_x*sqrt(2.*M_PI));

    pnl_vect_free(&Ci);
    pnl_vect_free(&log_Ci);
    pnl_vect_free(&log_A_func_T_ti);
    pnl_vect_free(&B_func_a_T_ti);
    pnl_vect_free(&B_func_b_T_ti);

    return Nominal*omega*P0_opmat*integral;

}


/*Payer Swaption payer_receiver=1, Receiver payer_receiver=-1*/
double cf_swaption_hw2d(ZCMarketData* ZCMarket, int payer_receiver, double Nominal, double periodicity, double option_maturity, double contract_maturity, double swaption_strike, double a, double b, double sigma, double eta, double rho)
{
    double xt, yt, price;
    xt=0.0;
    yt=0.0;

    // Transform the Hull/White model parameters to G2++ parameters.
    HW2dparams_to_G2dparams(a, b, &sigma, &eta, &rho);

    price = SWAPTION_g2(ZCMarket, payer_receiver, option_maturity, periodicity, contract_maturity, Nominal, swaption_strike, xt, yt, a, b, sigma, eta, rho);

    return price;
}

// Transform HW2d parameters into G2++ parameters
void HW2dparams_to_G2dparams(double a, double b, double *sigma, double *eta, double *rho)
{
    double sigma4, sigma3;

    sigma4 = (*eta) /fabs(a-b);
    sigma3 = sqrt(SQR((*sigma)) + SQR((*eta))/SQR(a-b) - 2*(*rho)*(*sigma)*(*eta)/fabs(a-b));

    (*rho) = ((*sigma)*(*rho) - sigma4) / sigma3;
    (*eta) = sigma4;
    (*sigma) = sigma3;
}

// Transform G2++ parameters into HW2d parameters
void G2dparams_to_HW2dparams(double a, double b, double *sigma, double *eta, double *rho)
{
    double sigma1, sigma2;

    sigma1 = sqrt(SQR((*sigma)) + SQR((*eta)) + 2*(*rho)*(*sigma)*(*eta));
    sigma2 = (*eta) * fabs(a-b);

    (*rho) = ((*sigma)*(*rho) + (*eta)) / sigma1;
    (*sigma) = sigma1;
    (*eta) = sigma2;
}


void HW2dparams_to_G2dparams_vect(PnlVect *HW2dparams, PnlVect *G2dparams)
{
    double a, sigma1, b, sigma2, rho;

    pnl_vect_resize(HW2dparams, 5);
    pnl_vect_resize(G2dparams, 5);

    a = GET(HW2dparams, 0);
    sigma1 = GET(HW2dparams, 1);
    b = GET(HW2dparams, 2);
    sigma2 = GET(HW2dparams, 3);
    rho = GET(HW2dparams, 4);
    HW2dparams_to_G2dparams(a, b, &sigma1, &sigma2, &rho);
    LET(G2dparams, 0) = a;
    LET(G2dparams, 1) = sigma1;
    LET(G2dparams, 2) = b;
    LET(G2dparams, 3) = sigma2;
    LET(G2dparams, 4) = rho;
}


