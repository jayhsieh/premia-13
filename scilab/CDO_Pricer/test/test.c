#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "mt19937.h"
#include <math.h>
#include "../src/cdo.h"

CDO			*laurent_gregory(void)
{
    double		x[2] = { 0., 5. };
//    double		a[] = { 0. };
    double		c[2] = {0.01, 0.01 };
    const double	delta = 0.4;
    const int		n_dates = 20;
    const double	t[20] = {0.25, 0.5, 0.75, 1., //
			       1.25, 1.5, 1.75, 2., //
			       2.25, 2.5, 2.75, 3., //
			       3.25, 3.5, 3.75, 4.,//
			       4.25, 4.5, 4.75, 5. };
    const int		n_tranches = 5;
    const double	tr[] = {0., 0.03, 0.06, 0.1, 1.};
    const int		n_comp = 100;
    company		**Co;
    double		*N;
    double		sN;
    int	i;

    sN = 0;
    N = malloc(n_comp * sizeof(double));
    for (i = 0; i < n_comp; i++){   
//	N[i] = 1./100.+0.005*gaussian();
	N[i] = 1./100.;
	sN += N[i];
    }
    Co = malloc(n_comp * sizeof(company*));
//    c[0] = 0.01;
    for (i = 0; i < n_comp; i++){   
	Co[i] = init_company_cov_cst(N[i]/sN, 1, x, c, delta);
//	Co[i] = init_company_cov_gauss(N[i]/sN, 1, x, c, delta, 0.2);
//	Co[i] = init_company_cov_unif(1./100., 1, x, c, 0.39, 0.41);
    }
    
    return (init_CDO(n_comp, Co, n_dates, t, n_tranches, tr));
}

int			main(void)
{
    CDO			*cdo;
    CDO			*hcdo;
    double		rho = 0.3;
    double		theta = 0.1964;
    copula		*cop;
    int			n_x = 100;
    grid		*x;
    grid		*t;
    int			j;
    int			jt;
    int			jx;
    int			jn;
    clock_t		start;
    clock_t		stop;
    cond_prob		*cp;
    double		**losses;
    double		**nd;
    grid		**meanloss;
    grid		**meanloss2;
    double		*dl;
    double		*pl;
    double		*hpl;
    double		*hdl;
    double		mcdl;
    double		mcpl;
    double		var_dl;
    double		var_pl;
    char		nom[20];
    FILE		*fich;
    int			n_mc;
    int			ntr;
    double		x2;
    double		xtaux[] = {0., 5., 10.};
    double		ytaux[] = {0.03, 0.15, 4.};
    step_fun		*rates, *r;
   /* 
    sf = init_cont_step_fun(3, sx, f_x);
    printf("%g\n", sf->f(sf, 3.)); 
*/
    start = clock();
    init_genrand(time(NULL));
    cdo = laurent_gregory();
    cop = init_gaussian_copula(rho); 
//    cop = init_clayton_copula(rho); 
//   cop = init_doublet_copula(rho, 5, 5); 
//    cop = init_nig_copula(rho, 1.2558, -0.2231); 
    x = init_hom_grid(MINDOUBLE, 0.6, 0.6/((double) n_x));
    rates = init_cont_linear_sf(2, xtaux, ytaux);
    for (x2 = 0.031; x2 < 3.8; x2+=0.001) 
	printf("%g\t%g\n", x2, inverse_sf(rates, x2));
    stop = clock();
    printf("%lu\t%lu\t%lu\n",start, stop, CLOCKS_PER_SEC);
    printf("%g\n", (double)(stop-start)/CLOCKS_PER_SEC);

 //   r = integrate_sf(rates);

/*    t = init_hom_grid(0.0125, 5., 0.0125);
    t = cdo->dates;
//    print_grid(cdo->dates);
    cp = init_cond_prob(cdo, cop, t);
    losses = lg_losses(cdo, cop, t, x, cp);
    for (jt = 0; jt < cdo->dates->size; jt++) {
        sprintf(nom, "lg_losses_%.2f_%.2f", rho, cdo->dates->data[jt]);
        fich = fopen(nom, "w+");
        for (jx = 0; jx < x->size; jx++)
            fprintf(fich, "%g\t%g\n", x->data[jx], losses[jt][jx]);
        fclose(fich);
    } 
    meanloss = mean_losses(cdo, t, x, losses);
    pl = payment_leg(cdo, rates, t, meanloss);
    dl = default_leg(cdo, rates, t, meanloss);
    ntr = cdo->n_tranches-1;
    for (j = 0; j < ntr; j++) {
	printf("[%.3f %.3f]\tPL: %.6f\tDL: %.6f\t Prix: %.2f\n", cdo->tr[j], cdo->tr[j+1], pl[j], dl[j], dl[j]/pl[j]/0.0001);
    }
    free(dl);
    free(pl);
    free(losses);
    free(cp);
    printf("\n");
*/
/*
ntr = cdo->n_tranches-1;
n_mc = 50000;
    printf("Payment Leg:\n");
    pl = mc_payment_leg(cdo, cop, rates, n_mc);
    printf("%f\t%f\n", pl[0], pl[4]);
    printf("%f\t%f\n", pl[1], pl[5]);
    printf("%f\t%f\n", pl[2], pl[6]);
    printf("%f\t%f\n", pl[3], pl[7]);
    printf("Default Leg:\n");
    dl = mc_default_leg(cdo, cop, rates, n_mc);
    printf("%f\t%f\n", dl[0], dl[4]);
    printf("%f\t%f\n", dl[1], dl[5]);
    printf("%f\t%f\n", dl[2], dl[6]);
    printf("%f\t%f\n", dl[3], dl[7]);
    for (j = 0; j < ntr; j++) {
	printf("[%.3f %.3f]\tPL: %.6f\tDL: %.6f\t Prix: %.2f\t[%.2f - %.2f]\n", cdo->tr[j], cdo->tr[j+1], pl[j], dl[j], dl[j]/pl[j]/0.0001,//
	(dl[j]-2.*sqrt(dl[j+ntr]/(double)n_mc))/(pl[j]+2.*sqrt(pl[j+ntr]/(double)n_mc))/0.0001,// 
	(dl[j]+2.*sqrt(dl[j+ntr]/(double)n_mc))/(pl[j]-2.*sqrt(pl[j+ntr]/(double)n_mc))/0.0001);
    }
    printf("Variable de controle...\n");
    hcdo = homogenize_CDO(cdo);
    t = cdo->dates;
//    t = init_hom_grid(0.025, 5., 0.025);
    cp = init_cond_prob(hcdo, cop, t);
    start = clock();
    losses = hw_losses_nh(hcdo, cop, t, x, cp);
    for (jt = 0; jt < cdo->dates->size; jt++) {
        sprintf(nom, "hw_losses_%.2f", cdo->dates->data[jt]);
        fich = fopen(nom, "w+");
        for (jx = 0; jx < x->size; jx++)
            fprintf(fich, "%g\t%g\n", x->data[jx], losses[jt][jx]);
        fclose(fich);
    }
    meanloss = mean_losses(hcdo, t, x, losses);
    hpl = payment_leg(hcdo, rates, t, meanloss);
    hdl = default_leg(hcdo, rates, t, meanloss);
    for (j = 0; j < ntr; j++) {
	printf("[%.3f %.3f]\tPL: %.6f\tDL: %.6f\t Prix: %.2f\n", cdo->tr[j], cdo->tr[j+1], hpl[j], hdl[j], hdl[j]/hpl[j]/0.0001);
    }
    free(meanloss);
    free(hpl);
    free(hdl);
    free(cp);
stop = clock();
printf("%f\n", (double)(stop - start) / CLOCKS_PER_SEC);

start = clock();
    cp = init_cond_prob(hcdo, cop, t);
    nd = hw_numdef(hcdo, cop, t, cp);
    for (jt = 0; jt < cdo->dates->size; jt++) {
        sprintf(nom, "hw_nd_%.2f", cdo->dates->data[jt]);
        fich = fopen(nom, "w+");
	for (jn = 0; jn < cdo->n_comp+1; jn++) {
            fprintf(fich, "%d\t%g\n", jn, nd[jt][jn]);
	}
        fclose(fich);
    }
    meanloss2 = mean_losses_from_numdef(hcdo, t, nd);
    hpl = payment_leg(hcdo, rates, t, meanloss2);
    hdl = default_leg(hcdo, rates, t, meanloss2);
    for (j = 0; j < ntr; j++) {
	printf("[%.3f %.3f]\tPL: %.6f\tDL: %.6f\t Prix: %.2f\n", cdo->tr[j], cdo->tr[j+1], hpl[j], hdl[j], hdl[j]/hpl[j]/0.0001);
    }
stop = clock();
printf("%f\n", (double)(stop - start) / CLOCKS_PER_SEC);
    
n_mc = 100000;
    printf("Payment Leg:\n");
    pl = mc_payment_vc_leg(cdo, cop, rates, n_mc);
    printf("%f\t%f\n", pl[0], pl[4]);
    printf("%f\t%f\n", pl[1], pl[5]);
    printf("%f\t%f\n", pl[2], pl[6]);
    printf("%f\t%f\n", pl[3], pl[7]);
    printf("Default Leg:\n");
    dl = mc_default_vc_leg(cdo, cop, rates, n_mc);
    printf("%f\t%f\n", dl[0], dl[4]);
    printf("%f\t%f\n", dl[1], dl[5]);
    printf("%f\t%f\n", dl[2], dl[6]);
    printf("%f\t%f\n", dl[3], dl[7]);
    for (j = 0; j < ntr; j++) {	
	printf("[%.3f %.3f]\tPL: %.6f\tDL: %.6f\t Prix: %.2f\t[%.2f - %.2f]\n", cdo->tr[j], cdo->tr[j+1], hpl[j]+pl[j], hdl[j]+dl[j], (hdl[j]+dl[j])/(hpl[j]+pl[j])/0.0001,//
	(hdl[j]+dl[j]-2.*sqrt(dl[j+ntr]/(double)n_mc))/(hpl[j]+pl[j]+2.*sqrt(pl[j+ntr]/(double)n_mc))/0.0001,// 
	(hdl[j]+dl[j]+2.*sqrt(dl[j+ntr]/(double)n_mc))/(hpl[j]+pl[j]-2.*sqrt(pl[j+ntr]/(double)n_mc))/0.0001);

    }
*/
/*    for (j = 0; j < ntr; j++) {
	printf("[%.3f %.3f]\tPL: %.6f\tDL: %.6f\t Prix: %.2f\t[%.2f - %.2f]\n", cdo->tr[j], cdo->tr[j+1], pl[j], dl[j], dl[j]/pl[j]/0.0001,//
	(dl[j]-2.*sqrt(dl[j+ntr]/(double)n_mc))/(pl[j]+2.*sqrt(pl[j+ntr]/(double)n_mc))/0.0001,// 
	(dl[j]+2.*sqrt(dl[j+ntr]/(double)n_mc))/(pl[j]-2.*sqrt(pl[j+ntr]/(double)n_mc))/0.0001);
    }
*/

/*    printf("%g\t%g\n", var_dl, var_pl);
    printf("[%.3f %.3f]\tPL: %.6f\tDL: %.6f\t Prix: %.2f\t[%.2f - %.2f]\n", cdo->A, cdo->B, mcpl, mcdl, mcdl/mcpl/0.0001,//
*/
/*
    cp = init_cond_prob(cdo, cop, t, v);
    losses = hw_losses_h(cdo, cop, t, x, v, cp);
    for (jt = 0; jt < cdo->dates->size; jt++) {
        sprintf(nom, "hw_losses_%.2f_%.2f", rho, cdo->dates->data[jt]);
        fich = fopen(nom, "w+");
        for (jx = 0; jx < x->size; jx++)
            fprintf(fich, "%g\t%g\n", x->data[jx], losses[jt][jx]);
        fclose(fich);
    }
    for (j = 0; j < cdo->n_tranches-1; j++) {
	cdo->A = cdo->tr[j];
	cdo->B = cdo->tr[j+1];
	meanloss = mean_losses(cdo, x, t, losses);
	pl = payment_leg(cdo, B, t, meanloss);
	dl = default_leg(B, t, meanloss);
	printf("[%.3f %.3f]\tPL: %.4f\tDL: %.4f\t Prix: %.2f\n", cdo->A, cdo->B, pl, dl, dl/pl/0.0001);
	free(meanloss);
    }
    free(losses);

    n_mc = 100000;
    cdo->A = cdo->tr[1];
    cdo->B = cdo->tr[2];
    dl = mc_default_leg(cdo, cop, B, n_mc);
    pl = mc_payment_leg(cdo, cop, B, n_mc);
    printf("[%.3f %.3f]\tPL: %.4f\tDL: %.4f\t Prix: %.2f\n", cdo->A, cdo->B, pl, dl, dl/pl/0.0001);
*/
/*    
    losses = mc_losses(cdo, cop, t, n_mc);
    for (j = 0; j < cdo->n_tranches-1; j++) {
	cdo->A = cdo->tr[j];
	cdo->B = cdo->tr[j+1];
	meanloss = mc_mean_losses(cdo, x, t, n_mc, losses);
	pl = payment_leg(cdo, B, t, meanloss);
	dl = default_leg(B, t, meanloss);
	printf("[%.3f %.3f]\tPL: %.4f\tDL: %.4f\t Prix: %.2f\n", cdo->A, cdo->B, pl, dl, dl/pl/0.0001);
	free(meanloss);
    }
  */
/*
  cp = init_cond_prob(cdo, cop, t, v);
    losses = lg_losses(cdo, cop, t, x, v, cp);
    for (jt = 0; jt < cdo->dates->size; jt++) {
        sprintf(nom, "hw_losses_%.2f_%.2f", rho, cdo->dates->data[jt]);
        fich = fopen(nom, "w+");
        for (jx = 0; jx < x->size; jx++)
            fprintf(fich, "%g\t%g\n", x->data[jx], losses[jt][jx]);
        fclose(fich);
    }

    for (j = 0; j < cdo->n_tranches-1; j++) {
	cdo->A = cdo->tr[j];
	cdo->B = cdo->tr[j+1];
	meanloss = mean_losses(cdo, x, t, losses);
	pl = payment_leg(cdo, B, t, meanloss);
	dl = default_leg(B, t, meanloss);
	printf("[%.2f %.2f]\tPL: %.4f\tDL: %.4f\t Prix: %.2f\n", cdo->A, cdo->B, pl, dl, dl/pl/0.0001);
	free(meanloss);
    }
*/
    return (0);
}

/*printf("temps: %g\n", ((double) (end - start)) / CLOCKS_PER_SEC);*/
/*    free(losses);*/
/*   cp = init_cond_prob(cdo, cop, t, v);
    printf("essai\n");
    nd = hw_numdef(cdo, cop, t, v, cp);
    for (jt = 0; jt < cdo->dates->size; jt++) {
        sprintf(nom, "hw_nd_%f", cdo->dates->data[jt]);
        fich = fopen(nom, "w+");
	for (jn = 0; jn < cdo->n_comp+1; jn++) {
            fprintf(fich, "%d\t%g\n", jn, nd[jt][jn]);
	}
        fclose(fich);
    }
    return (0);
}
*/
 /*
start = clock();
    lg_ws = lg_create_workspace(cdo, cop, t, x, v);
    losses = lg_losses(cdo, cop, t, x, v, lg_ws); */
/*    for (jt = 0; jt < cdo->dates->size; jt++) {*/
/*        sprintf(nom, "lg_losses_%f", cdo->dates->data[jt]);*/
/*        fich = fopen(nom, "w+");*/
/*        for (jx = 0; jx < x->size; jx++)*/
/*            fprintf(fich, "%g\t%g\n", x->data[jx], losses[jt][jx]);*/
/*        fclose(fich);*/
/*    }*/
/*
    cdo->A = 0;
    cdo->B = 0.03;
    meanloss = mean_losses(cdo, x, t, losses);
    pl = payment_leg(cdo, B, t, meanloss);
    dl = default_leg(B, t, meanloss);
    printf("[%.2f %.2f]\t %f\t %f\t %f\n", cdo->A, cdo->B, pl, dl, dl/pl/0.0001);
    cdo->A = 0.03;
    cdo->B = 0.06;
    meanloss = mean_losses(cdo, x, t, losses);
    pl = payment_leg(cdo, B, t, meanloss);
    dl = default_leg(B, t, meanloss);
    printf("[%.2f %.2f]\t %f\t %f\t %f\n", cdo->A, cdo->B, pl, dl, dl/pl/0.0001);
    cdo->A = 0.06;
    cdo->B = 0.10;
    meanloss = mean_losses(cdo, x, t, losses);
    pl = payment_leg(cdo, B, t, meanloss);
    dl = default_leg(B, t, meanloss);
    printf("[%.2f %.2f]\t %f\t %f\t %f\n", cdo->A, cdo->B, pl, dl, dl/pl/0.0001);
    cdo->A = 0.10;
    cdo->B = 1.;
    meanloss = mean_losses(cdo, x, t, losses);
    pl = payment_leg(cdo, B, t, meanloss);
    dl = default_leg(B, t, meanloss);
    printf("[%.2f %.2f]\t %f\t %f\t %f\n", cdo->A, cdo->B, pl, dl, dl/pl/0.0001);
end = clock();
*/
/*printf("temps: %.g\n", ((double) (end - start)) / CLOCKS_PER_SEC);*/
  //  free(losses);

/*    cdo->A = 0;*/
/*    cdo->B = 0.03;*/
/*    pl = payment_leg(cdo, B, x, t, losses);*/
/*    dl = default_leg2(cdo, B, x, t, losses);*/
/*    printf("PL: %g\tDL: %g\t Prix: %g\n", pl, dl, dl/pl/0.0001);*/
/*    cdo->A = 0.1;*/
/*    cdo->B = 1.;*/
/*    pl = payment_leg(cdo, B, x, t, losses);*/
/*    dl = default_leg2(cdo, B, x, t, losses);*/
/*    printf("PL: %g\tDL: %g\t Prix: %g\n", pl, dl, dl/pl/0.0001);*/

/*    cjl = lg_cond_jump_losses(cdo, cop, t, x, v, ws);*/
/*    printf("OK \n");*/
/*    cdo->A = 0;*/
/*    cdo->B = 0.03;*/
/*    pl = payment_leg(cdo, B, x, t, losses);*/
/*    pl += lg_accrued_mPL(cdo, cop, t, x, v, ws, cjl, B);*/
/*    dl = lg_default_leg(cdo, cop, t, x, v, ws, cjl, B);*/
/*    printf("PL: %g\tDL: %g\t Prix: %g\n", pl, dl, dl/pl/0.0001);*/
/*    cdo->A = 0.03;*/
/*    cdo->B = 0.1;*/
/*    pl = payment_leg(cdo, B, x, t, losses);*/
/*    pl += lg_accrued_mPL(cdo, cop, t, x, v, ws, cjl, B);*/
/*    dl = lg_default_leg(cdo, cop, t, x, v, ws, cjl, B);*/
/*    printf("PL: %g\tDL: %g\t Prix: %g\n", pl, dl, dl/pl/0.0001);*/
/*    cdo->A = 0.1;*/
/*    cdo->B = 1.;*/
/*    pl = payment_leg(cdo, B, x, t, losses);*/
/*    pl += lg_accrued_mPL(cdo, cop, t, x, v, ws, cjl, B);*/
/*    dl = lg_default_leg(cdo, cop, t, x, v, ws, cjl, B);*/
/*    printf("PL: %g\tDL: %g\t Prix: %g\n", pl, dl, dl/pl/0.0001);*/
/*
    free_grid(t);
    free_grid(v);
    free_grid(x);
    free(cdo);
    free(cop);

    return (0);
}

*/
