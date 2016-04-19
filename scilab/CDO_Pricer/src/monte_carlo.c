#include        "mt19937.h"
#include        "cdo.h"

static int      compute_default(const CDO       *cdo,
                                copula          *cop,
                                int         *ind,
                                double          *tau)
{
  double      tau_jn;
  int         n_def;
  int         jn;
  int         jk;
    
  cop->generate(cop);
  n_def = 0;
  for (jn = 0; jn < cdo->n_comp; jn++) {
    if (cop->compute_default_time(cop, cdo->C[jn]->H, &tau_jn)) {
      jk = n_def-1;
      while ((jk >= 0) && (tau_jn < tau[jk])) {
        ind[jk+1] = ind[jk];
        tau[jk+1] = tau[jk];
        jk--;
      }
      ind[jk+1] = jn;
      tau[jk+1] = tau_jn;
      n_def++;
    }
  }

  return (n_def);
}

double          *mc_default_one_leg(const CDO       *cdo,
                                    copula      *cop,
                                    const step_fun  *rates,
                                    int         *ind,
                                    double      *tau)
{
  double      *dl;
  double      *phi_losses;
  double      losses;
  double      act;
  double      new_phi_losses;
  int         n_def;
  int         jk;
  int         jtr;
  int         jt;
    
  n_def = compute_default(cdo, cop, ind, tau);
  dl = malloc((cdo->n_tranches-1) * sizeof(double));
  phi_losses = malloc((cdo->n_tranches-1) * sizeof(double));
  for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) { 
    dl[jtr] = 0;
    phi_losses[jtr] = 0.;
  }
  losses = 0;
  jt = 0;
  for (jk = 0; jk < n_def; jk++) {
    losses += cdo->C[ind[jk]]->nominal * (1. - RECOVERY(ind[jk]));
    act = exp(- compute_sf(rates, tau[jk]));    
    for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) {
      new_phi_losses = pp(losses, cdo->tr[jtr]) - pp(losses, cdo->tr[jtr+1]);
      dl[jtr] += act * (new_phi_losses - phi_losses[jtr]);
      phi_losses[jtr] = new_phi_losses;
    }
  }
  free(phi_losses);

  return (dl); 
}

double          *mc_payment_one_leg(const CDO       *cdo,
                                    copula      *cop,
                                    const step_fun  *rates,
                                    int         *ind,
                                    double      *tau)
{
  double      losses;
  double      *pl;
  double      act;
  double      new_phi_losses;
  double      *phi_losses;
  double      t;
  int         jtr;
  int         n_def;
  int         jk;
  int         jt;
    
  n_def = compute_default(cdo, cop, ind, tau);
  pl = malloc((cdo->n_tranches-1) * sizeof(double));
  phi_losses = malloc((cdo->n_tranches-1) * sizeof(double));
  for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) { 
    pl[jtr] = 0;
    phi_losses[jtr] = 0.;
  }
  losses = 0;
  jk = 0;
  t = 0;
  for (jt = 0; jt < cdo->dates->size; jt++) {
    while ( (tau[jk] >= t)&&(tau[jk] < cdo->dates->data[jt])&&(jk<n_def) ) {
      losses += cdo->C[ind[jk]]->nominal * (1. - RECOVERY(ind[jk]));
      jk++;
    }
    act = exp(- compute_sf(rates, cdo->dates->data[jt])); 
    for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) {
      new_phi_losses = pp(losses, cdo->tr[jtr]) - pp(losses, cdo->tr[jtr+1]);
      pl[jtr] += act * (cdo->tr[jtr+1] - cdo->tr[jtr] - new_phi_losses) * (cdo->dates->data[jt] - t); 
    }
    t = cdo->dates->data[jt];
  }
  losses = 0;
  jt = 0;
  for (jk = 0; jk < n_def; jk++) {
    while (tau[jk] > cdo->dates->data[jt]) jt++;
    t = (jt == 0) ? 0. : cdo->dates->data[jt-1];
    losses += cdo->C[ind[jk]]->nominal * (1. - RECOVERY(ind[jk]));
    act = exp(- compute_sf(rates, tau[jk]));    
    for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) {
      new_phi_losses = pp(losses, cdo->tr[jtr]) - pp(losses, cdo->tr[jtr+1]);
      pl[jtr] += act * (new_phi_losses - phi_losses[jtr]) * (tau[jk] - t);
      phi_losses[jtr] = new_phi_losses;
    }
  }
  free(phi_losses);

  return (pl); 
}
    
double          *mc_default_vc_one_leg(const CDO        *cdo,
                                       copula      *cop,
                                       const step_fun  *rates,
                                       int         *ind,
                                       double      *tau)
{
  double      *dl;
  double      *phi_losses;
  double      *phi_losses_vc;
  double      losses;
  double      losses_vc;
  double      act;
  double      new_phi_losses;
  double      new_phi_losses_vc;
  double      nominal;
  double      delta;
  int         n_def;
  int         jk;
  int         jtr;
  int         jt;
  int	      jc;

  n_def = compute_default(cdo, cop, ind, tau);
  dl = malloc((cdo->n_tranches-1) * sizeof(double));
  phi_losses = malloc((cdo->n_tranches-1) * sizeof(double));
  phi_losses_vc = malloc((cdo->n_tranches-1) * sizeof(double));
  for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) { 
    dl[jtr] = 0;
    phi_losses[jtr] = 0.;
    phi_losses_vc[jtr] = 0.;
  }
 
  nominal = 0;
  delta = 0;
  for (jc = 0; jc < cdo->n_comp; jc++){
    nominal += cdo->C[jc]->nominal;
    delta += cdo->C[jc]->mean_delta;
  }
  nominal /= (double) cdo->n_comp;
  delta /= (double) cdo->n_comp;

  losses = 0;
  losses_vc = 0;
  jt = 0;
  for (jk = 0; jk < n_def; jk++) {
    losses += cdo->C[ind[jk]]->nominal * (1. - RECOVERY(ind[jk]));
    losses_vc += nominal * (1. - delta);
    act = exp(- compute_sf(rates, tau[jk]));    
    for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) {
      new_phi_losses = pp(losses, cdo->tr[jtr]) - pp(losses, cdo->tr[jtr+1]);
      new_phi_losses_vc = pp(losses_vc, cdo->tr[jtr]) - pp(losses_vc, cdo->tr[jtr+1]);
      dl[jtr] += act * (new_phi_losses - phi_losses[jtr] - (new_phi_losses_vc - phi_losses_vc[jtr]));
      phi_losses[jtr] = new_phi_losses;
      phi_losses_vc[jtr] = new_phi_losses_vc;
    }
  }
  free(phi_losses);
  free(phi_losses_vc);

  return (dl); 
}
    
double          *mc_payment_vc_one_leg(const CDO        *cdo,
                                       copula      *cop,
                                       const step_fun  *rates,
                                       int         *ind,
                                       double      *tau)
{
  double      losses;
  double      losses_vc;
  double      *pl;
  double      act;
  double      new_phi_losses;
  double      new_phi_losses_vc;
  double      *phi_losses;
  double      *phi_losses_vc;
  double      t;
  double      nominal;
  double      delta;
  int         jtr;
  int         n_def;
  int         jk;
  int         jt;
  int	      jc;
    
  n_def = compute_default(cdo, cop, ind, tau);
  pl = malloc((cdo->n_tranches-1) * sizeof(double));
  phi_losses = malloc((cdo->n_tranches-1) * sizeof(double));
  phi_losses_vc = malloc((cdo->n_tranches-1) * sizeof(double));
  for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) { 
    pl[jtr] = 0;
    phi_losses[jtr] = 0.;
    phi_losses_vc[jtr] = 0.;
  }
  
  nominal = 0;
  delta = 0;
  for (jc = 0; jc < cdo->n_comp; jc++){
    nominal += cdo->C[jc]->nominal;
    delta += cdo->C[jc]->mean_delta;
  }
  nominal /= (double) cdo->n_comp;
  delta /= (double) cdo->n_comp;

  losses = 0;
  losses_vc = 0;
  jk = 0;
  t = 0;
  for (jt = 0; jt < cdo->dates->size; jt++) {
    while ( (tau[jk] >= t)&&(tau[jk] < cdo->dates->data[jt])&&(jk<n_def) ) {
      losses += cdo->C[ind[jk]]->nominal * (1. - RECOVERY(ind[jk]));
      losses_vc += nominal * (1. - delta); 
      jk++;
    }
    act = exp(- compute_sf(rates, cdo->dates->data[jt])); 
    for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) {
      new_phi_losses = pp(losses, cdo->tr[jtr]) - pp(losses, cdo->tr[jtr+1]);
      new_phi_losses_vc = pp(losses_vc, cdo->tr[jtr]) - pp(losses_vc, cdo->tr[jtr+1]);
      pl[jtr] += act * (new_phi_losses_vc - new_phi_losses) * (cdo->dates->data[jt] - t); 
    }
    t = cdo->dates->data[jt];
  }
  losses = 0;
  losses_vc = 0;
  jt = 0;
  for (jk = 0; jk < n_def; jk++) {
    while (tau[jk] > cdo->dates->data[jt]) jt++;
    t = (jt == 0) ? 0. : cdo->dates->data[jt-1];
    losses += cdo->C[ind[jk]]->nominal * (1. - RECOVERY(ind[jk]));
    losses_vc += nominal * (1. - delta); 
    act = exp(- compute_sf(rates, tau[jk]));    
    for (jtr = 0; jtr < cdo->n_tranches-1; jtr++) {
      new_phi_losses = pp(losses, cdo->tr[jtr]) - pp(losses, cdo->tr[jtr+1]);
      new_phi_losses_vc = pp(losses_vc, cdo->tr[jtr]) - pp(losses_vc, cdo->tr[jtr+1]);
      pl[jtr] += act * (new_phi_losses - phi_losses[jtr] - (new_phi_losses_vc - phi_losses_vc[jtr])) * (tau[jk] - t);
      phi_losses[jtr] = new_phi_losses;
      phi_losses_vc[jtr] = new_phi_losses_vc;
    }
  }
  free(phi_losses);
  free(phi_losses_vc);

  return (pl); 
}
    
double          *mc_generic_leg(const CDO       *cdo,
                                copula          *cop,
                                const step_fun      *rates,
                                const int       n_mc,
                                mc_one_leg      *one_leg)
{
  double      *leg;
  double      **stock;
  double      *tau;
  int         *ind;
  int	      jnc;
  int         jmc;
  int         jtr;
  int         ntr = cdo->n_tranches-1;

  leg = malloc(2 * (ntr) * sizeof(double));
  stock = malloc(n_mc * sizeof(double*));
  tau = malloc(cdo->n_comp * sizeof(double));
  ind = malloc(cdo->n_comp * sizeof(int));
  for (jtr = 0; jtr < 2*ntr; jtr++)
    leg[jtr] = 0;
  for (jnc = 0; jnc < cdo->n_comp; jnc++) {
    tau[jnc] = 0;
    ind[jnc] = 0;
  }
  for (jmc = 0; jmc < n_mc; jmc++) {
    stock[jmc] = one_leg(cdo, cop, rates, ind, tau);
    for (jtr = 0; jtr < ntr; jtr++)
      leg[jtr] += stock[jmc][jtr];
  }
  for (jtr = 0; jtr < ntr; jtr++)
    leg[jtr] /= (double) n_mc;
  for (jmc = 0; jmc < n_mc; jmc++) {
    for (jtr = 0; jtr < ntr; jtr++)
      leg[ntr+jtr] += (stock[jmc][jtr] - leg[jtr]) * (stock[jmc][jtr] - leg[jtr]);
  }
  for (jtr = 0; jtr < ntr; jtr++)
    leg[ntr+jtr] /= ((double) n_mc - 1.);
  free(ind);
  free(tau);
  free(stock);

  return (leg);
}

double          *mc_default_leg(const CDO       *cdo,
                                copula          *cop,
                                const step_fun      *rates,
                                const int       n_mc)
{
  return (mc_generic_leg(cdo, cop, rates, n_mc, mc_default_one_leg));
}

double          *mc_payment_leg(const CDO       *cdo,
                                copula          *cop,
                                const step_fun      *rates,
                                const int       n_mc)
{
  return (mc_generic_leg(cdo, cop, rates, n_mc, mc_payment_one_leg));
}

double          *mc_default_vc_leg(const CDO        *cdo,
                                   copula          *cop,
                                   const step_fun      *rates,
                                   const int       n_mc)
{
  return (mc_generic_leg(cdo, cop, rates, n_mc, mc_default_vc_one_leg));
}

double          *mc_payment_vc_leg(const CDO        *cdo,
                                   copula          *cop,
                                   const step_fun      *rates,
                                   const int       n_mc)
{
  return (mc_generic_leg(cdo, cop, rates, n_mc, mc_payment_vc_one_leg));
}


