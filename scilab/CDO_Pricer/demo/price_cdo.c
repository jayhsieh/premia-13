#include "../src/cdo.h"

void            print_losses(const CDO      *cdo,
                             const grid     *t,
                             const grid     *x,
                             double* const  *losses)
{
  int         jn;
  int         jt;
  int         jtdates;
  char        nom[50];
  FILE        *fich;
    
  for (jt = 0, jtdates = 0; jt < t->size; jt++) {
    if (t->data[jt] == cdo->dates->data[jtdates]) {
      sprintf(nom, ".losses_%i", jtdates);
      fich = fopen(nom, "w+");
      for (jn = 0; jn < x->size; jn++) {
        fprintf(fich, "%g\t%g\n", x->data[jn], losses[jt][jn]);
      }
      fclose(fich);
      jtdates++;
    }
  }
}

void            print_numdef(const CDO      *cdo,
                             const grid     *t,
                             double* const  *numdef)
{
  int         jn;
  int         jt;
  int         jtdates;
  char        nom[50];
  FILE        *fich;
    
  for (jt = 0, jtdates = 0; jt < t->size; jt++) {
    if (t->data[jt] == cdo->dates->data[jtdates]) {
      sprintf(nom, ".numdef_%i", jtdates);
      fich = fopen(nom, "w+");
      for (jn = 0; jn < cdo->n_comp+1; jn++) {
        fprintf(fich, "%i\t%g\n", jn, numdef[jt][jn]);
      }
      fclose(fich);
      jtdates++;
    }
  }
}

extern int      price_cdo(const int    *n_comp, 
                          const double *nominal,
                          const int     n_dates,
                          const double *dates,
                          const int     n_tranches,
                          const double *tr,
                          const double *intensity,
                          const int     n_rates,
                          const double *x_rates,
                          const double *y_rates,
                          const int    *t_recovery,
                          const double *recovery,
                          const int    *t_copula,
                          const double *p_copula,
                          const int    *t_method,
                          const int    *p_method,
                          double       *price,
                          double       *def_leg,
                          double       *pay_leg)
{
  int n=50;
  prod* produit;
  produit=malloc(sizeof(prod));
  produit->nominal=nominal[0];
  produit->nb=n_comp[0];
  produit->recov=recovery[0];
  produit->maturite=n_dates*dates[0];
  double *p;
  double *rho;
  double *tab;
  int e;
  double     **tab_sadd;
  double     ***U;

  int         n_sub = 1;
  double      x_intensity[2];
  double      y_intensity[2]; 
  double      i[n_sub];
  int         jt;
  int         jtr;
  int         n_mc;
  int         size_price;
  step_fun    *rates;
  company     **Co;
  CDO         *cdo;
  CDO         *hcdo;
  copula      *cop;
  grid        *x = NULL;
  grid        *t = NULL;
  cond_prob   *cp = NULL;
  double      **numdef = NULL;
  double      **losses = NULL;
  grid        **meanloss = NULL;
  int         jn;
  double      *dl;
  double      *pl;
  double      *hdl;
  double      *hpl;
  double      sum_nominal;
    
  Co = malloc(*n_comp * sizeof(company*));
  x_intensity[0] = 0.;
  x_intensity[1] = dates[n_dates-1];
  switch (*t_recovery) {
  case 1 :
    for (jn = 0; jn < *n_comp; jn++) { 
      y_intensity[0] = intensity[jn];
      y_intensity[1] = intensity[jn];
      Co[jn] = init_company_cov_cst(nominal[jn], n_sub, x_intensity, y_intensity, recovery[0]);
    }
    break;
  case 2 :
    for (jn = 0; jn < *n_comp; jn++) { 
      y_intensity[0] = intensity[jn];
      y_intensity[1] = intensity[jn];
      Co[jn] = init_company_cov_unif(nominal[jn], n_sub, x_intensity, y_intensity, recovery[0], recovery[1]);
    }
    break;
  case 3 :
    for (jn = 0; jn < *n_comp; jn++) { 
      y_intensity[0] = intensity[jn];
      y_intensity[1] = intensity[jn];
      Co[jn] = init_company_cov_gauss(nominal[jn], n_sub, x_intensity, y_intensity, recovery[0], recovery[1]);
    }
    break;
  }
  cdo = init_CDO(*n_comp, Co, n_dates, dates, n_tranches, tr);
  switch (*t_copula) {
  case 1 :
    cop = init_gaussian_copula(p_copula[0]);
    break;
  case 2 :
    cop = init_clayton_copula(p_copula[0]);
    break;
  case 3 :
    cop = init_nig_copula(p_copula[0], p_copula[1], p_copula[2]);
    break;
  case 4:
    cop=  init_student_copula( p_copula[0],p_copula[1]);
    break;
  case 5:
    cop=  init_double_t_copula( p_copula[0],p_copula[1],p_copula[2]);
    break;
  }
  rates = init_cont_linear_sf(n_rates-1, x_rates, y_rates);
 
  produit->rate=compute_sf(rates,1);

  switch (*t_method)
    {
    case 1 :
      if(*t_copula!=4){
        t = init_fine_grid(cdo->dates, p_method[0]);
        cp = init_cond_prob(cdo, cop, t);
        numdef = hw_numdef(cdo, cop, t, cp);
        print_numdef(cdo, t, numdef);
        meanloss = mean_losses_from_numdef(cdo, t, numdef);
        pl = payment_leg(cdo, rates, t, meanloss);
        dl = default_leg(cdo, rates, t, meanloss);
      }
      else{
        t = init_fine_grid(cdo->dates, p_method[0]);
        cp = init_cond_prob(cdo, cop, t);
        numdef = hw_numdef1(cdo, cop, t, cp);
        print_numdef(cdo, t, numdef);
        meanloss = mean_losses_from_numdef(cdo, t, numdef);
        pl = payment_leg(cdo, rates, t, meanloss);
        dl = default_leg(cdo, rates, t, meanloss);
      }
      break;
    case 2 :
      if(*t_copula!=4){
        t = init_fine_grid(cdo->dates, p_method[0]);
        cp = init_cond_prob(cdo, cop, t);
        numdef = lg_numdef(cdo, cop, t, cp);
        print_numdef(cdo, t, numdef);
        meanloss = mean_losses_from_numdef(cdo, t, numdef);
        pl = payment_leg(cdo, rates, t, meanloss);
        dl = default_leg(cdo, rates, t, meanloss);
      }
      else  {
        t = init_fine_grid(cdo->dates, p_method[0]);
        cp = init_cond_prob(cdo, cop, t);
        numdef = lg_numdef1(cdo, cop, t, cp);
   
        print_numdef(cdo, t, numdef);
        meanloss = mean_losses_from_numdef(cdo, t, numdef);
        pl = payment_leg(cdo, rates, t, meanloss);
        dl = default_leg(cdo, rates, t, meanloss);
      }
      break;
    case 3 :
      if(*t_copula!=4) {
        t = init_fine_grid(cdo->dates, p_method[0]);
        x = init_hom_grid(MINDOUBLE, (1.-recovery[0]), (1.-recovery[0])/(double) p_method[1]);
        cp = init_cond_prob(cdo, cop, t);
        losses = hw_losses_h(cdo, cop, t, x, cp);
        print_losses(cdo, t, x, losses);
        meanloss = mean_losses(cdo, t, x, losses);
        pl = payment_leg(cdo, rates, t, meanloss);
        dl = default_leg(cdo, rates, t, meanloss);
      }
      else{
        t = init_fine_grid(cdo->dates, p_method[0]);
        x = init_hom_grid(MINDOUBLE, (1.-recovery[0]), (1.-recovery[0])/(double) p_method[1]);
        cp = init_cond_prob(cdo, cop, t);
        losses = hw_losses_h1(cdo, cop, t, x, cp);
        print_losses(cdo, t, x, losses);
        meanloss = mean_losses(cdo, t, x, losses);
        pl = payment_leg(cdo, rates, t, meanloss);
        dl = default_leg(cdo, rates, t, meanloss);
      }
     
      break;
    case 4 :
      if(*t_copula!=4){
        t = init_fine_grid(cdo->dates, p_method[0]);
        x = init_hom_grid(MINDOUBLE, (1.-recovery[0]), (1.-recovery[0])/(double) p_method[1]);
        cp = init_cond_prob(cdo, cop, t);
        losses = hw_losses_nh(cdo, cop, t, x, cp);
        print_losses(cdo, t, x, losses);
        meanloss = mean_losses(cdo, t, x, losses);
        pl = payment_leg(cdo, rates, t, meanloss);
        dl = default_leg(cdo, rates, t, meanloss);
      }
      else{  
        t = init_fine_grid(cdo->dates, p_method[0]);
        x = init_hom_grid(MINDOUBLE, (1.-recovery[0]), (1.-recovery[0])/(double) p_method[1]);
        cp = init_cond_prob(cdo, cop, t);
        losses = hw_losses_nh1(cdo, cop, t, x, cp);
        print_losses(cdo, t, x, losses);
        meanloss = mean_losses(cdo, t, x, losses);
        pl = payment_leg(cdo, rates, t, meanloss);
        dl = default_leg(cdo, rates, t, meanloss);
      }

      break;
    case 5 : 
      if(*t_copula!=4){   
        t = init_fine_grid(cdo->dates, p_method[0]);
        x = init_hom_grid(MINDOUBLE, (1.-recovery[0]), (1.-recovery[0])/(double) p_method[1]);
        cp = init_cond_prob(cdo, cop, t);
        losses = lg_losses(cdo, cop, t, x, cp);
        print_losses(cdo, t, x, losses);
        meanloss = mean_losses(cdo, t, x, losses);
        pl = payment_leg(cdo, rates, t, meanloss);
        dl = default_leg(cdo, rates, t, meanloss);
      }
      else{
        t = init_fine_grid(cdo->dates, p_method[0]);
        x = init_hom_grid(MINDOUBLE, (1.-recovery[0]), (1.-recovery[0])/(double) p_method[1]);
        cp = init_cond_prob(cdo, cop, t);
        losses = lg_losses1(cdo, cop, t, x, cp);
        print_losses(cdo, t, x, losses);
        meanloss = mean_losses(cdo, t, x, losses);
        pl = payment_leg(cdo, rates, t, meanloss);
        dl = default_leg(cdo, rates, t, meanloss);
      }
      break;
    case 6 :
      n_mc = p_method[0]; 
      pl = mc_payment_leg(cdo, cop, rates, n_mc);
      dl = mc_default_leg(cdo, cop, rates, n_mc);
   
      break;
    case 7 :
      if(*t_copula!=4) {
        n_mc = p_method[0];
        t = init_fine_grid(cdo->dates, p_method[1]);
        hcdo = homogenize_CDO(cdo);
        cp = init_cond_prob(hcdo, cop, t);
        numdef = lg_numdef(hcdo, cop, t, cp);
        meanloss = mean_losses_from_numdef(hcdo, t, numdef);
        hpl = payment_leg(hcdo, rates, t, meanloss);
        hdl = default_leg(hcdo, rates, t, meanloss);
        pl = mc_payment_vc_leg(cdo, cop, rates, n_mc);
        dl = mc_default_vc_leg(cdo, cop, rates, n_mc);
        for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) {
          pl[jtr] = hpl[jtr] + pl[jtr];
          dl[jtr] = hdl[jtr] + dl[jtr];
        }
      }
      else
        {
          n_mc = p_method[0];
          t = init_fine_grid(cdo->dates, p_method[1]);
          hcdo = homogenize_CDO(cdo);
          cp = init_cond_prob(hcdo, cop, t);
          numdef = lg_numdef1(hcdo, cop, t, cp);
          meanloss = mean_losses_from_numdef(hcdo, t, numdef);
          hpl = payment_leg(hcdo, rates, t, meanloss);
          hdl = default_leg(hcdo, rates, t, meanloss);
          pl = mc_payment_vc_leg(cdo, cop, rates, n_mc);
          dl = mc_default_vc_leg(cdo, cop, rates, n_mc);
          for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) {
            pl[jtr] = hpl[jtr] + pl[jtr];
            dl[jtr] = hdl[jtr] + dl[jtr];
          }
        }
      break;
    
    case 8:
    
      if((*t_copula!=4)&&(*t_copula!=2))
        {
          pl=malloc((cdo->n_tranches-1)*sizeof(double));
          dl=malloc((cdo->n_tranches-1)*sizeof(double));
          t = init_fine_grid(cdo->dates, p_method[0]);
          cp = init_cond_prob(cdo, cop, t);
          tab_sadd=malloc((cdo->n_tranches-1)*sizeof(double*));
          U=malloc((cdo->n_tranches-1)*sizeof(double**));
 
          for(jtr=0;jtr<cdo->n_tranches-1;jtr++){
            tab_sadd[jtr]=malloc((t->size)*sizeof(double));
            U[jtr]=malloc((t->size)*sizeof(double*)); 
          }
          for(jtr=0;jtr<cdo->n_tranches-1;jtr++){
            for(jt=0;jt<t->size;jt++){
              U[jtr][jt]=malloc((cop->size)*sizeof(double**));
            }
          }
  
          U=Uoptimal(cdo,cop,t,cp); 
          tab_sadd=saddlepoint(cdo,cop,t,cp,U);
  
          dl=default_leg_sadd(cdo,rates,t,tab_sadd,cop);      
          pl=payment_leg_sadd(cdo,rates,t,tab_sadd,cop);
        }
      else
        {
          printf("NON TREATED CASE\n");
          return(0);
        }
 
      break;


    case 9:

      rho=malloc(5*sizeof(double));
      tab=malloc(5*sizeof(double));
   
      for(e=0;e<5;e++){
        tab[e]=p_method[e];
      }

      pl=malloc((4.0+cdo->n_tranches)*sizeof(double));
      dl=malloc((4.0+cdo->n_tranches)*sizeof(double)); 
  
      rho=Base_correl(produit,tab,p_method[5],p_copula[0]);

      for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) {
        produit->att=tr[jtr];
        produit->det=tr[jtr+1];
      
        pl[jtr] =pay_leg_base(produit,tab,p_method[5],p_copula[0]);
        dl[jtr] =dl_leg_base(produit,tab,p_method[5],p_copula[0]);
        
      }
      for(jtr=cdo->n_tranches-1;jtr<cdo->n_tranches+4;jtr++){
        pl[jtr]=rho[jtr-cdo->n_tranches+1];
      }
      if(p_copula[0]==1){
        for(jtr=0;jtr<4;jtr++){
          dl[jtr+cdo->n_tranches-1]=0.015+0.03*jtr;
        }
        dl[cdo->n_tranches+3]=0.17;
      }
      else if(p_copula[0]==2){
        dl[cdo->n_tranches-1]=0.015;
        dl[cdo->n_tranches]=0.05;      
        dl[cdo->n_tranches+1]=0.065;
        dl[cdo->n_tranches+2]=0.1125;
        dl[cdo->n_tranches+3]=0.225;
      }

      free(tab);
      free(rho);      
     
      
      break;

    case 10:
  
      p=malloc(n*sizeof(double));
      tab=malloc(5*sizeof(double));

      for(e=0;e<5;e++){
        tab[e]=p_method[e];
      }
      p=probaimpl(produit,tab,p_method[5],n,p_copula[0]);
      pl=malloc((n+cdo->n_tranches-1)*sizeof(double));
      dl=malloc((n+cdo->n_tranches-1)*sizeof(double));
   

      for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) {
        produit->att=tr[jtr];
        produit->det=tr[jtr+1];
        pl[jtr] =pay_leg_impl(produit,p,n) ;
        dl[jtr] =dl_leg_impl(produit,p,n)  ;
      }
      dl[cdo->n_tranches-1]=p[0];
      pl[cdo->n_tranches-1]=0;   
      for(jtr=cdo->n_tranches;jtr<cdo->n_tranches+n-1;jtr++){
        pl[jtr]=(1+jtr-cdo->n_tranches)*0.1*1./(n-1);
        dl[jtr]=p[1+jtr-cdo->n_tranches]*(n-1)/0.1;
      }


      free(tab);
      free(p); 

      break;

    }

  
    
  size_price = n_tranches-1;
  if ((*t_method==6 )) size_price = 2*(n_tranches-1);
  if ((*t_method==7 )) size_price = 2*(n_tranches-1);
  if((*t_method ==9 )) size_price = (n_tranches+4.0);
  if((*t_method ==10 )) size_price = (n_tranches+49);
  
  for (jtr = 0; jtr < size_price; jtr++) {
    pay_leg[jtr] = (double) pl[jtr];
    def_leg[jtr] = (double) dl[jtr];
    price[jtr] = (def_leg[jtr] / (pay_leg[jtr])) * 10000.;
  }

  
  free_cdo(cdo);

  if (numdef != NULL) {
    for (jt = 0; jt < t->size; jt++) 
      free(numdef[jt]);
    free(numdef);
  }
  if (losses != NULL) {
    for (jt = 0; jt < t->size; jt++) 
      free(losses[jt]);
    free(losses);
  }
  if (cp != NULL) free_cond_prob(cp);
  if (meanloss != NULL) {
    for (jtr = 0; jtr < n_tranches-1; jtr++)
      free_grid(meanloss[jtr]);
    free(meanloss);
  }
  if (t != NULL) free_grid(t);
  if (x != NULL) free_grid(x);
  free(pl);
  free(dl);
  
  return (0);
}
