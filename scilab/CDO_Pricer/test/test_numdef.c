#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include "../src/cdo.h"

int			test_numdef(const int	n_comp, 
				    const double	*N,
				    const double	*intensity_cst,
				    const double	recovery,
				    const int	n_dates,
				    const double	*dates,
				    const int	n_tranches,
				    const double	*tr,
				    const double	rho,
				    const int	n_v)
{
    int	n_sub = 1;
    double		x_int[] = { 0. };
    int	jn;
    company		**Co;
    CDO			*cdo;
    copula		*gc;
    grid		*v;
    grid		*x;
    grid		*t;
    cond_prob		*cp;
    double		**nd;
    

    Co = malloc(n_comp * sizeof(company*));
    for (jn = 0; jn < n_comp; jn++)
	Co[jn] = init_company_cov_cst(N[jn], n_sub, x_int, intensity_cst, recovery);
    init_CDO(n_comp, Co, n_dates, dates, n_tranches, tr);
    gc = init_gaussian_copula(rho);
    v = init_hom_grid(-12., 12., 24./((double) n_v));
//    x = init_hom_grid(MINDOUBLE, 0.6, 0.6/((double) n_x));
    t = init_hom_grid(0.05, 5., 0.05);
    cp = init_cond_prob(cdo, gc, t, v);
    nd = hw_numdef(cdo, gc, t, v, cp);
    /*    
    for (jt = 0; jt < cdo->dates->size; jt++) {
        sprintf(nom, "hw_nd_%f", cdo->dates->data[jt]);
        fich = fopen(nom, "w+");
	for (jn = 0; jn < cdo->n_comp+1; jn++) {
            fprintf(fich, "%d\t%g\n", jn, nd[jt][jn]);
	}
        fclose(fich);
    }
    */

    return (0);
}

