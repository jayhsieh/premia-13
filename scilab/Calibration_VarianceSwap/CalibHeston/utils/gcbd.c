#include "optim.h"

static  int
optim_zgcbd (opt_simul simul, int *n, double *binf, double *bsup, double *x,
	     double *f, double *g, double *zero, int *napmax, int *itmax,
	     int *indgc, int *ibloc, int *nfac, int *imp, int *io,
	     double *epsx, double *epsf, double *epsg, double *dir,
	     double *df0, double *diag, double *x2, 
	     opt_simul_data *optim_data,
	     double *y, double *s, double *z__, double *ys,
	     double *zs, int *nt, int *index, double *wk1, double *wk2,
	     double *alg, int *ialg, char *nomf, long int nomf_len);

static int optim_gcp (int *n, int *index, int *indic, int *np, int *nt, double *y,
		      double *s, double *z__, double *ys, double *zs, double *diag,
		      double *b, double *x, double *d__, double *g, double *eps);

static int optim_calbx (int *n, int *index, int *indic, int *nt, int *np, double *y,
			double *s, double *ys, double *z__, double *zs, double *x,
			double *diag, double *bx);

static double optim_rednor (int *n, double *binf, double *bsup, double *x,
			    double *epsx, double *g);

static int optim_shanph (double *diag, int *n, int *nt, int *np, double *y, double *s,
			 double *ys, double *scal, int *index, int *io, int *imp);

static int optim_relvar (int *ind, int *n, double *x, double *binf, double *bsup,
			 double *x2, double *g, double *diag, int *imp, int *io,
			 int *ibloc, int *izag, int *iter, int *nfac, int *irit);

static int optim_bfgsd (double *diag, int *n, int *nt, int *np, double *y, double *s,
			double *ys, double *condm, double *param, double *zero,
			int *index);

static int  optim_majz (int *n, int *np, int *nt, double *y, double *s, double *z__,
			double *ys, double *zs, double *diag, int *index);

/*!but 
 *    algorithme de minimisation d une fonction reguliere sous 
 *    contraintes de borne 
 *!methode 
 *    methode de bfgs a memoire limitee + projection 
 *!origine 
 *    f. bonnans inria juin 1985 
 *    Copyright INRIA 
 * 
 *!sous programmes (modulopt) 
 *    proj rlbd majysa majz calbx gcp relvar bfgsd shanph 
 *!liste d' appel 
 *    indgc   indicateur de gcbd                                  es 
 *      en entree =1 standard 
 *                =2 dh et indic initialises au debut de trav et itrav 
 *                   ifac,f,g initialises 
 *      en sortie 
 *       si < 0 incapacite de calculer un point meilleur que le point initial 
 *       si = 0 arret demande par l utilisateur 
 *       si > 0 on fournit un point meilleur que le point de depart 
 *       = -14 insuffisance memoire 
 *       = -13 indgc non egal a zero ou 1 en entree 
 *       = -12 zero,epsg ou df0 non strict. positifs 
 *       = -11 n,napmax,itmax ou io non strict. positifs 
 *       < -10 parametres d entree non convenables 
 *       = -6  arret lors du calcul de la direction de descente et iter=1 
 *       = -5  arret lors du calcul de l approximation du hessien  iter=1 
 *       = -3  anomalie de simul : indic negatif en un point ou 
 *             f et g ont ete precedemment calcules 
 *       = -2  echec de la recherche lineaire a la premiere iteration 
 *       = -1  f non definie au point initial 
 *       =  1  arret sur epsg 
 *       =  2            epsf 
 *       =  3            epsx 
 *       =  4            napmax 
 *       =  5            itmax 
 *       =  6  pente dans la direction opposee au gradient trop petite 
 *       =  7  arret lors du calcul de la direction de descente 
 *       =  8  arret lors du calcul de l approximation du hessien 
 *       = 10  arret par echec de la recherche lineaire , cause non precisee 
 *       = 11  idem avec indsim < 0 
 *       = 12            un pas trop petit proche d un pas trop grand 
 *                       ceci peut resulter d une erreur dans le gradient 
 *       = 13            trop grand nombre d appels dans une recherche lineaire 
 * 
 *    simul  subroutine permettant de calculer f et g (norme modulopt) 
 *    n dim de x                                                 e 
 *    x variables a optimiser (controle)                          es 
 *    f valeur du critere                                         s 
 *    g gradient de f                                             s 
 *    imp si =0 pas d impression 
 *            1  impressions en debut etfin dexecution 
 *            2  3 lignes a chaque iteration 
 *            >=3 nombreuses impressions    e 
 *    io numero fichier sortie            e 
 *    zero  proche zero machine                                             e 
 *    napmax nombre maximum d appels de simul                               e 
 *    itmax nombre maximum d iterations de gcbd                             e 
 *    epsf critere arret sur f            e 
 *    epsg arret si sup a norm2(g+)/n     e 
 *    epsx vect dim n precision sur x     e 
 *    df0>0 decroissance f prevue         e 
 *    binf,bsup bornes inf,sup,de dim n                          e 
 *    nfac nombre de variables non bloquees a l optimum          s 
 *    vect,ivect vecteurs de travail de dim nvect,nivect 
 *    izs,rzs,dzs : cf normes modulopt         es 
 * 
 *! 
 *        signification de quelques variables internes 
 * 
 *    {y}={g1}-{g0}                                        l (locale) 
 *    {s}={x1}-{x0}                                        l 
 *    {z}=[b]*{s}. [b] est une estimation de hessien       l 
 *    ys=<y>*{s}                                           l 
 *    zs=<z>*{s}                                           l 
 *    diag approximation diagonale du hessien  es 
 *    si indgc=0 diag initialise a ******************* 
 *    puis remis a jour par bfgs diagonal 
 *    nt: le nombre de deplacements pris en compte          l 
 *    index(nt) repere les vect y,s,z                       l 
 *    wk1,wk2: vecteurs de travail de dim n                 l 
 *    ibloc vect dim n  ; si x(i) est bloque, ibloc(i)=iteration de blocage ; 
 *    si x(i) est libre, ibloc(i)=-1*(iteration de deblocage) 
 *    irit: irit=1, si relachement de vars a l'iter courante, 0 sinon 
 *    ired: ired=1 decision de redemarrage, 0 sinon 
 *    alg(1)=param 
 *    alg(2)=condmax 
 *    alg(6)=eps 
 *    alg(9)=tetaq ( critere de redemarrage) 
 *    ialg(1)       correction de powell sur y si (y,s)trop petit 
 *         0:       sans correction de powell 
 *         1:       avec correction 
 *    ialg(2)       mise a jour de diag par bfgsd 
 *         0:       pas de remise a jour 
 *         1:       remise a jour de diag par bfgsd 
 *    ialg(3)       mise a l'echelle par methode de shanno-phua 
 *         0:       pas de mise a l'echelle 
 *         1:       mise a l'echelle a toutes les iterations 
 *         2:       mise a l'echelle a la 2ieme iteration seulement 
 *    ialg(4):      memorisation pour choix iterations 
 *         0:                sans memorisation 
 *         1:      avec memorisation 
 *    ialg(5):      memorisation par variable 
 *         0:      sans memorisation 
 *         1:       avec memorisation 
 *    ialg(6):      choix des iterations de relachement 
 *         1:      relachement a toutes les iterations 
 *         2:      relachement si decroissance g norme gradient 
 *         10:     relachement si decroissance critere % iter.initcycle 
 *         11:      relachement si decroissance critere % decroissance cycle 
 *    ialg(7):      choix des variables a relacher 
 *         1:       methode de bertsekas modifiee 
 *    ialg(8):      choix de la direction de descente 
 *         3:      qn sans memoire: nt derniers deplacements 
 *         4:      redemarrage sans accumulation 
 *         5:      redemarrage avec accumulation 
 *    ialg(9):     critere de redemarrage 
 *         2:       redemarrage si fact. ou defact. 
 *         10:     decroissance critere % decroissance iter_init du cycle 
 *         11:     decroissance critere % decroissance totale du cycle 
 *         12:      diminution de znglib d un facteur alg(9)=tetaq 
 *    eps0 sert a partitionner les variables 
 *    ceps0 utilise dans le calcul de eps0 
 *    izag nombre d iterations min pendant lesquelles une var reste bloquee 
 *    nap nombre d appels de simul 
 *    iter iteration courante 
 *    ind indicateur de simul 
 *    icv memoire entree indgc 
 *    np  param utilise pour la gestion de vect. nb de vect courant. 
 *    lb  param utilise pour la gestion de vect. 1er place libre. 
 *    nb  param utilise pour la gestion de vect. 
 *       nb=2 si l'algorithme utilise est qn sans memoire +redem +pas acc 
 *       nb=1 sinon 
 * 
 * 
 *    initialisation des parametres 
 */


int optim_gcbd (int *indgc, opt_simul simul, char *nomf, int *n, double *x, double *f,
		double *g, int *imp, int *io, double *zero, int *napmax,
		int *itmax, double *epsf, double *epsg, double *epsx, double *df0,
		double *binf, double *bsup, int *nfac, double *vect, int *nvect,
		int *ivect, int *nivect,
		opt_simul_data *optim_data,
		long int nomf_len)
{
  int i__1;
  double d__1, d__2;

  /* Local variables */
  int ialg[15], nfin, ndir, i__, ndiag;
  double aa;
  int ii, nd, ng, ns, nt, nindic, ny, nz, nindex, nx2;
  double alg[15];
  int nys, nzs;

  /* Parameter adjustments */
  --bsup;
  --binf;
  --epsx;
  --g;
  --x;
  --vect;
  --ivect;


  /* Function Body */
  nt = 2;
  alg[0] = 1e-5;
  alg[1] = 1e6;
  alg[5] = .5;
  alg[8] = .5;
  /* 
   */
  ialg[0] = 1;
  ialg[1] = 0;
  ialg[2] = 2;
  ialg[3] = 0;
  ialg[4] = 0;
  ialg[5] = 2;
  ialg[6] = 1;
  ialg[7] = 4;
  ialg[8] = 12;
  /* 
   *    verification des entrees 
   *Computing MIN 
   */
  i__1 = Min (*n, *napmax);
  ii = Min (i__1, *itmax);
  if (ii > 0)
    {
      goto L10;
    }
  *indgc = -11;
  if (*imp > 0)
    {
      Sciprintf("gcbd: return with indg=%d\n",*indgc);
    }
  return 0;
 L10:
  /*Computing MIN 
   */
  d__1 = Min (*zero, *epsg);
  aa = Min (d__1, *df0);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L11: */
      /*Computing MIN 
       */
      d__1 = aa, d__2 = epsx[i__];
      aa = Min (d__1, d__2);
    }
  if (aa > 0.)
    {
      goto L12;
    }
  *indgc = -12;
  if (*imp > 0)
    { 
      Sciprintf("gcbd: return with indg=%d\n",*indgc);
    }
  return 0;
 L12:
  /* 
   *    decoupage de la memoire 
   */
  ny = 1;
  ns = nt * *n + ny;
  nz = nt * *n + ns;
  nys = nt * *n + nz;
  nzs = nt + nys;
  nd = nt + nzs;
  ng = *n + nd;
  nx2 = *n + ng;
  ndir = *n + nx2;
  ndiag = *n + ndir;
  nfin = *n + ndiag;
  /* 
   */
  if (nfin > *nvect)
    {   
      Sciprintf("gcbd: insufficient memory, nvect=%d and should be %d\n",nfin,*nvect);
      *indgc = -14;
      return 0;
    }
  /* 
   */
  nindic = 1;
  nindex = *n + nindic;
  nfin = nt + nindex;
  if (nfin > *nivect)
    {
      Sciprintf("gcbd: insufficient memory, nivect=%d and should be %d\n",nfin,*nivect);
      *indgc = -14;
      return 0;
    }
  /* 
   */
  optim_zgcbd ( simul, n, &binf[1], &bsup[1], &x[1], f, &g[1], zero,
		napmax, itmax, indgc, &ivect[nindic], nfac, imp, io, &epsx[1],
		epsf, epsg, &vect[ndir], df0, &vect[ndiag], &vect[nx2],
		optim_data, &vect[ny], &vect[ns], &vect[nz],
		&vect[nys], &vect[nzs], &nt, &ivect[nindex], &vect[nd],
		&vect[ng], alg, ialg, nomf, 6L);
  return 0;
}





static  int
optim_zgcbd (opt_simul simul, int *n, double *binf, double *bsup, double *x,
	     double *f, double *g, double *zero, int *napmax, int *itmax,
	     int *indgc, int *ibloc, int *nfac, int *imp, int *io,
	     double *epsx, double *epsf, double *epsg, double *dir,
	     double *df0, double *diag, double *x2, 
	     opt_simul_data *optim_data,
	     double *y, double *s, double *z__, double *ys,
	     double *zs, int *nt, int *index, double *wk1, double *wk2,
	     double *alg, int *ialg, char *nomf, long int nomf_len)
{
  int c__1 = 1;
  /* System generated locals */
  int y_dim1, y_offset, s_dim1, s_offset, z_dim1, z_offset, i__1;
  double d__1, d__2;

  /* Local variables */
  double diff, difg, scal;
  int ired;
  double diri;
  int nred, izag;
  int napm;
  double teta;
  int iter, irit;
  double tmax, ceps0;
  int izag1, napm1;
  double teta1, znog0;
  int i__;
  double t, condm, param;
  int icycl, napav, indrl;
  double tetaq;
  double epsxi, tproj;
  int indgc1;
  double dfred1, param1, dfrit1, aa;
  int lb, nb;
  double fn, difred;
  int np;
  double xi, sy, epsgcp;
  int indsim;
  double znglib, difrit;
  double zngred;
  int iresul;
  double zngrit, ys1, amd, amf;
  int ind;
  double dfp;
  int nap, ifp, irl, inp;
  double bss, zng, zrl;
  int imp1;
  double eps0, bss2;



  /*    Copyright INRIA 
   * 
   * 
   */
  /* Parameter adjustments */
  --wk2;
  --wk1;
  --x2;
  --diag;
  --dir;
  --epsx;
  --ibloc;
  --g;
  --x;
  --bsup;
  --binf;
  --index;
  --zs;
  --ys;
  z_dim1 = *nt;
  z_offset = z_dim1 + 1;
  z__ -= z_offset;
  s_dim1 = *nt;
  s_offset = s_dim1 + 1;
  s -= s_offset;
  y_dim1 = *nt;
  y_offset = y_dim1 + 1;
  y -= y_offset;
  --alg;
  --ialg;

  /* Function Body */
  if (*imp >= 4)
    {
      Sciprintf("dans gcbd. algorithme utilise: \n");
      if (ialg[1] == 1)
	{
	  Sciprintf("\temploi correction de powell \n");
	}
      if (ialg[2] == 1)
	{
	  Sciprintf("\tmise a jour de diag par la methode bfgs\n");
	}
      if (ialg[3] == 1)
	{
	  Sciprintf("\tmise a echelle de diag par methode de shanno-phua\n");
	}
      if (ialg[3] == 2)
	{
	  Sciprintf("\tmise a echelle de diag seulement a la 2e iter\n");
	}
      if (ialg[4] == 1)
	{
	  Sciprintf("\tmemorisation pour choix iteration \n");
	}
      if (ialg[5] == 1)
	{
	  Sciprintf("\tmemorisation par variable\n");
	}
      if (ialg[6] == 1)
	{
	  Sciprintf("\trelachememt de variables a toutes les iteration\n");
	}
      if (ialg[6] == 2)
	{
	  Sciprintf("\trelachement de vars si decroissance g_norme\n");
	}
      if (ialg[6] == 10)
	{
	  Sciprintf("\trelachement de vars si dec f %% iter_init du cycle\n");
	}
      if (ialg[6] == 11)
	{
	  Sciprintf("\trelachement de vars si dec f %% dec du cycle\n");
	}
      if (ialg[7] == 1)
	{
	  Sciprintf("\tchoix de vars a relacher par bertsekas modifiee\n");
	}
      if (ialg[8] == 1)
	{
	  Sciprintf("\tchoix de dir descente par methode de gradient\n");
	}
      if (ialg[8] == 2)
	{
	  Sciprintf("\tchoix de dir descente par methode qn\n");
	}
      if (ialg[8] == 3)
	{
	  Sciprintf("\tchoix de dir descente par qn sans memoire.nt depl\n");
	}
      if (ialg[8] == 4)
	{
	  Sciprintf("\tchoix de dir descente par qn -mem,redem,sans acc.\n");
	}
      if (ialg[8] == 5)
	{
	  Sciprintf("\tchoix de dir descente par qn -mem,redem,avec acc.\n");
	}
      if (ialg[9] == 2)
	{
	  Sciprintf( "\tredem si relachement de vars\n");
	}
      if (ialg[9] == 10)
	{
	  Sciprintf("\tredem si dec f %% dec iter_init du cycle\n");
	}
      if (ialg[9] == 11)
	{
	  Sciprintf("\tredem si dec f %% dec totale du cycle.\n");
	}
      if (ialg[9] == 12)
	{
	  Sciprintf("\tredem si diminution du gradient des var libres d un,facteur,%11.4f\n",
		    alg[9]);
	}
    }
  /* 
   *    section 1  initialisations 
   *    irl nombre de rech lin 'lentes' 
   *    nred nombre de redemarrage de la direction de descente 
   *    icycl nombre de cycles de minimisation 
   * 
   */
  epsgcp = 1e-5;
  indsim = 4;
  indrl = 1;
  irl = 0;
  irl = 0;
  nred = 1;
  icycl = 1;
  nap = 0;
  /* 
   */
  iresul = 1;
  optim_proj (n, &binf[1], &bsup[1], &x[1]);
  indsim = 4;
  (*simul) (&indsim, n, &x[1], f, &g[1],optim_data);
  ++nap;
  if (indsim > 0)
    {
      goto L99;
    }
  *indgc = -1;
  if (indsim == 0)
    {
      *indgc = 0;
    }
  if (*imp > 0)
    {
      Sciprintf( "gcbd : retour avec indgc=%d\n", (*indgc));
    }
  goto L900;
 L99:
  ceps0 = 20.;
  eps0 = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L100: */
      eps0 += epsx[i__];
    }
  eps0 = ceps0 * eps0 / *n;
  /* 
   *    calcul de zng 
   */
  znog0 = optim_rednor (n, &binf[1], &bsup[1], &x[1], &epsx[1], &g[1]);
  zng = znog0;
  zngrit = znog0;
  zngred = znog0;
  /* 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L130: */
      ibloc[i__] = 0;
    }
  izag = 3;
  izag1 = izag;
  nap = 0;
  iter = 0;
  scal = 1.;
  *nfac = *n;
  np = 0;
  lb = 1;
  nb = 2;
  if (ialg[8] == 3)
    {
      nb = 1;
    }
  i__1 = *nt;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L140: */
      index[i__] = i__;
    }
  tetaq = alg[9];
  condm = alg[2];
  param = alg[1];
  indgc1 = *indgc;
  /*    si indgc=0 on init diag a k*ident puis scal a it=2 
   * 
   */
  if (*indgc == 1 || *indgc >= 100)
    {
      goto L150;
    }
  if (*indgc == 2)
    {
      goto L180;
    }
  *indgc = -13;
  if (*imp > 0)
    {
      Sciprintf( "gcbd : retour avec indgc=%d\n", (*indgc));
    }
  goto L900;
  /* 
   */
 L150:
  /*    on initialise diag par approximation quadratique 
   *    df0 decroissance prevue . si mod quad df0=((dh)-1g,g)/2 
   *    et on cherche dh diag de la forme cst/(dx)**2 
   *    donc cst=som((g(i)*(dx))**2))/(2*df0) 
   */
  sy = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L160: */
      /*Computing 2nd power 
       */
      d__1 = g[i__] * epsx[i__];
      sy += d__1 * d__1;
    }
  sy /= *df0 * 2.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L170: */
      /*Computing 2nd power 
       */
      d__1 = epsx[i__];
      diag[i__] = (sy + *zero) / (d__1 * d__1 + *zero);
    }
 L180:
  /* 
   * 
   *    bouclage 
   */
 L200:
  ++iter;
  *indgc = 1;
  if (iter > *itmax)
    {
      *indgc = 5;
      goto L900;
    }
  /* L201: */
  if (*imp >= 2)
    {
      Sciprintf("(dans gcbd  iter=%d, f=%15.7f\n",iter,(*f));
    }
  if (iter == 1)
    {
      irit = 1;
      goto L301;
    }
  /* 
   */
  optim_majysa (n, nt, &np, &y[y_offset], &s[s_offset], &ys[1], &lb, &g[1],
		&x[1], &wk2[1], &wk1[1], &index[1], &ialg[1], &nb);
  inp = index[np];
  /* 
   * 
   *    correction powell sur y si (y,s) trop petit 
   */
  if (ialg[1] != 1)
    {
      goto L290;
    }
  param1 = 1. - param;
  bss = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L260: */
      /*Computing 2nd power 
       */
      d__1 = s[inp + i__ * s_dim1];
      bss += diag[i__] * (d__1 * d__1);
    }
  bss2 = param * bss;
  if (ys[inp] > bss2)
    {
      goto L290;
    }
  if (*imp > 2)
    {
      Sciprintf("gcbd. emploi correction powell (y,s)=%11.4f\n",ys[inp]);
    }
  teta = param1 * bss / (bss - ys[inp]);
  teta1 = 1. - teta;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L274: */
      y[inp + i__ * y_dim1] =
	teta * y[inp + i__ * y_dim1] + teta1 * diag[i__] * s[inp +
							     i__ * s_dim1];
    }
  ys[inp] = bss2;
  /*    verif correction powell (facultatif; faire go to 300) 
   */
  ys1 = C2F(ddot) (n, &s[inp + s_dim1], &c__1, &y[inp + y_dim1], &c__1);
  ys1 = (d__1 = bss2 - ys1, Abs (d__1)) / bss2;
  if (*imp > 2)
    {
      Sciprintf("erreur relative correction powell =%11.4f\n",ys1);
    }
  /* 
   *mise a jour de diag 
   */
 L290:
  if (ialg[2] == 1)
    {
      optim_bfgsd (&diag[1], n, nt, &np, &y[y_offset], &s[s_offset], &ys[1],
		   &condm, &param, zero, &index[1]);
    }
  /* 
   */
  if (ialg[3] == 1 || (ialg[3] == 2 && iter == 2))
    {
      optim_shanph (&diag[1], n, nt, &np, &y[y_offset], &s[s_offset], &ys[1],
		    &scal, &index[1], io, imp);
    }
  /* 
   */
  optim_majz (n, &np, nt, &y[y_offset], &s[s_offset], &z__[z_offset], &ys[1],
	      &zs[1], &diag[1], &index[1]);
  /* 
   *    section 3 determination des variables libres et bloquees 
   */
  /* L300: */
  /*    -----decision de relachement a l'iteration courante 
   *         relachement si irit=1 (sinon irit=0) 
   */
  irit = 0;
  if (ialg[6] == 1)
    {
      irit = 1;
    }
  if (ialg[6] == 2 && znglib <= alg[6] * zngrit)
    {
      irit = 1;
    }
  if (ialg[6] == 10 && diff <= dfrit1 * alg[6])
    {
      irit = 1;
    }
  if (ialg[6] == 11 && diff <= difrit * alg[6])
    {
      irit = 1;
    }
  if (irit == 1)
    {
      ++nred;
    }
  /*   ----choix des variables a relacher 
   */
  imp1 = *imp;
 L301:
  if (ialg[7] == 1)
    {
      optim_relvar (&ind, n, &x[1], &binf[1], &bsup[1], &x2[1], &g[1],
		    &diag[1], imp, io, &ibloc[1], &izag, &iter, nfac, &irit);
    }
  /* 
   * 
   *    section 4 expression de dir 
   */
  if (np == 0)
    {
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  dir[i__] = -g[i__] / diag[i__];
	  /* L400: */
	}
    }
  else
    {
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  dir[i__] = -scal * g[i__];
	  /* L410: */
	}
      optim_gcp (n, &index[1], &ibloc[1], &np, nt, &y[y_offset], &s[s_offset],
		 &z__[z_offset], &ys[1], &zs[1], &diag[1], &g[1], &dir[1],
		 &wk1[1], &wk2[1], &epsgcp);
    }
  /* 
   *    section 5  redemarrage 
   * 
   */
  if (ialg[8] == 4 || ialg[8] == 5)
    {
      ired = 0;
      if (ialg[9] == 2 && ind == 1)
	{
	  ired = 1;
	}
      if (ialg[9] == 10 && diff < dfred1 * tetaq)
	{
	  ired = 1;
	}
      if (ialg[9] == 11 && diff < difred * tetaq)
	{
	  ired = 1;
	}
      if (ialg[9] == 12 && znglib <= tetaq * zngred)
	{
	  ired = 1;
	}
      if (ired == 1)
	{
	  ++icycl;
	  np = 0;
	  lb = 1;
	  if (*imp > 2)
	    {
	      Sciprintf("gcbd  redemarrage. icycl=%d\n",	 icycl);
	    }
	}
    }
  /* 
   *    section 6 annulation de d(i) , i dans ib 
   */
  if (ialg[6] == 1)
    {
      goto L640;
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L630: */
      if (ibloc[i__] > 0)
	{
	  dir[i__] = 0.;
	}
    }
 L640:
  /* 
   *    recherche lineaire 
   *    conservation de x et g dans wk1 et wk2 
   */
  C2F(dcopy) (n, &x[1], &c__1, &wk1[1], &c__1);
  C2F(dcopy) (n, &g[1], &c__1, &wk2[1], &c__1);
  /*    calcul de la derivee dans la direction dir 
   */
  ifp = 0;
  fn = *f;
  znog0 = zng;
 L702:
  dfp = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      epsxi = epsx[i__];
      xi = x[i__];
      diri = dir[i__];
      if (xi - binf[i__] <= epsxi && diri < 0.)
	{
	  dir[i__] = 0.;
	}
      /* L710: */
      if (bsup[i__] - xi <= epsxi && diri > 0.)
	{
	  dir[i__] = 0.;
	}
    }
  dfp = C2F(ddot) (n, &g[1], &c__1, &dir[1], &c__1);
  if (-dfp > 0.)
    {
      goto L715;
    }
  if (ifp == 1)
    {
      *indgc = 6;
      goto L900;
    }
  /*    restauration dir 
   */
  if (*imp >= 3)
    {
      Sciprintf("gcbd : restauration dir ; fp=%11.4f,zero=%11.4f\n",	dfp,	(*zero));
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L712: */
      dir[i__] = -scal * g[i__];
    }
  ifp = 1;
  goto L702;
 L715:
  /*    pas initial suivant idee fletcher 
   */
  t = diff * -2. / dfp;
  if (iter == 1)
    {
      t = *df0 * -2. / dfp;
    }
  tmax = 1e10;
  t = Min (t, tmax);
  /*Computing MAX 
   */
  d__1 = t, d__2 = *zero * 1e10;
  t = Max (d__1, d__2);
  napm = 15;
  napm1 = nap + napm;
  if (napm1 > *napmax)
    {
      napm1 = *napmax;
    }
  napav = nap;
  amd = .7;
  amf = .1;
  /* 
   */
  optim_rlbd (&indrl, n,  simul, &x[1], &binf[1], &bsup[1], f, &dfp, &t,
	      &tmax, &dir[1], &g[1], &tproj, &amd, &amf, imp, io, zero, &nap,
	      &napm1, &x2[1],optim_data);
  if (*imp > 2)
    {
      Sciprintf("gcbd retour mlibd indrl=%d, pas=%11.4f, f=%11.4f\n",	indrl,t,(*f));
    }
  if (nap - napav >= 5)
    {
      ++irl;
    }
  if (indrl >= 10)
    {
      indsim = 4;
      ++nap;
      (*simul) (&indsim, n, &x[1], f, &g[1],optim_data);
      if (indsim <= 0)
	{
	  *indgc = -3;
	  if (indsim == 0)
	    {
	      *indgc = 0;
	    }
	  if (*imp > 0)
	    {
	      Sciprintf( "gcbd : retour avec indgc=%d\n", (*indgc));
	    }
	  goto L900;
	}
    }
  if (indrl <= 0)
    {
      *indgc = 10;
      if (indrl == 0)
	{
	  *indgc = 0;
	}
      if (indrl == -3)
	{
	  *indgc = 13;
	}
      if (indrl == -4)
	{
	  *indgc = 12;
	}
      if (indrl <= -1000)
	{
	  *indgc = 11;
	}
      if (*imp > 0)
	{

	  Sciprintf( "gcbd : retour avec indgc=%d\n",  (*indgc));
	}
      goto L900;
    }
  if (*imp >= 5)
    {
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  /* L760: */
	  if (*imp > 2)
	    {
	      Sciprintf( "gcbd i=%d, xgd ,%11.4f,%11.4f,%11.4f\n",  i__,
			 x[i__],
			 g[i__],
			 dir[i__]);
	    }
	}
    }
  /* 
   */
  if (nap < *napmax)
    {
      goto L758;
    }
  if (*imp > 0)
    {
      Sciprintf( "gcbd max appels simul\n");

    }
  *indgc = 4;
  goto L900;
 L758:
  /* 
   *    section 8 test de convergence 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if ((d__1 = x[i__] - wk1[i__], Abs (d__1)) > epsx[i__])
	{
	  goto L806;
	}
      /* L805: */
    }
  if (*imp > 0)
    {
      Sciprintf( "gcbd. retour apres convergence sur x\n");
    }
  *indgc = 3;
  goto L900;
  /*    calcul grad residuel,norme l2 
   */
 L806:
  difg = optim_rednor (n, &binf[1], &bsup[1], &x[1], &epsx[1], &g[1]);
  diff = fn - *f;
  if (*imp >= 2)
    {
      Sciprintf("gcbd. epsg=%11.4f, difg=%11.4f, epsf=%11.4f, diff=%11.4f, nap=%d\n",
		(*epsg),
		difg,
		(*epsf),
		diff,
		nap);
    }
  /* 
   */
  if (diff <= *epsf)
    {
      *indgc = 2;
      goto L900;
    }
  if (difg <= *epsg)
    {
      *indgc = 1;
      goto L900;
    }
  /* 
   *    -----mise a jour de difrit,dfrit1,difred,dfred1 
   */
  if (irit == 1)
    {
      difrit = diff;
      dfrit1 = diff;
    }
  else
    {
      difrit += diff;
    }
  if (ired == 1)
    {
      difred = diff;
      dfred1 = diff;
    }
  else
    {
      difred += diff;
    }
  /* 
   */
  znglib = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (ibloc[i__] > 0)
	{
	  goto L884;
	}
      aa = g[i__];
      if (x[i__] - binf[i__] <= epsx[i__])
	{
	  aa = Min (0., aa);
	}
      if (bsup[i__] - x[i__] <= epsx[i__])
	{
	  aa = Max (0., aa);
	}
      /*Computing 2nd power 
       */
      d__1 = aa;
      znglib += d__1 * d__1;
    L884:
      ;
    }
  znglib = sqrt (znglib);
  if (ired == 1)
    {
      zngred = znglib;
    }
  if (irit == 1)
    {
      zngrit = znglib;
    }
  goto L200;
  /* 
   *     fin des calculs 
   */
 L900:
  if (indrl == 0)
    {
      *indgc = 0;
    }
  if (*indgc == 1 && indrl <= 0)
    {
      *indgc = indrl;
    }
  if (*imp > 0)
    {
      Sciprintf( "gcbd : retour avec indgc=%d\n",	 (*indgc));
    }
  if (*imp >= 1 && (double) indrl <= *zero)
    {
      Sciprintf("arret impose par la recherche lineaire. cf notice rlbd,/, indicateur de rlbd=%d\n",
		indrl);
    }
  if (*imp >= 1)
    {
      Sciprintf("gcbd f,norme grad=%11.4f,difg=%11.4f ,nap=%d,iter=%d,indgc=%d\n",
		(*f),
		difg,
		nap,
		iter,
		(*indgc));
    }
  /* 
   *    autres impressions finales 
   */
  if (indgc1 < 100)
    {
      return 0;
    }
  zrl = 0.;
  if (iter > 0)
    {
      zrl = (double) nap / (double) iter;
    }
  /* L2000: */

  Sciprintf( "gcbd (%s,%11.4f,%11.4f,%d,%d,%6.2f,%d)\n",    nomf,
	     (*f),
	     difg,
	     nap,
	     iter,
	     zrl,
	     irl);
  return 0;
}	


/* 
 *    methode de gradient preconditionne appliquee a l'equation 
 *    [a]*{x}={b}. ici [a] est definie par les vecteurs ({y}(i), 
 *    {s}(i), {z}(i), i=1,np). 
 *    Copyright INRIA 
 * 
 * 
 *    initialisation 
 */

static int optim_gcp (int *n, int *index, int *indic, int *np, int *nt, double *y,
		      double *s, double *z__, double *ys, double *zs, double *diag,
		      double *b, double *x, double *d__, double *g, double *eps)
{
  /* System generated locals */
  int y_dim1, y_offset, s_dim1, s_offset, z_dim1, z_offset, i__1;

  /* Local variables */
  double beta;
  int iter, i__;
  int itmax;
  double s0, s1, s2, dg, ro, d2a, eps0, eps1;

  /* Parameter adjustments */
  --g;
  --d__;
  --x;
  --b;
  --diag;
  --indic;
  --zs;
  --ys;
  z_dim1 = *nt;
  z_offset = z_dim1 + 1;
  z__ -= z_offset;
  s_dim1 = *nt;
  s_offset = s_dim1 + 1;
  s -= s_offset;
  y_dim1 = *nt;
  y_offset = y_dim1 + 1;
  y -= y_offset;
  --index;

  /* Function Body */
  eps0 = 1e-5;
  eps1 = 1e-5;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L100;
	}
      x[i__] = -b[i__] / diag[i__];
    L100:
      ;
    }
  /* 
   */
  optim_calbx (n, &index[1], &indic[1], nt, np, &y[y_offset], &s[s_offset],
	       &ys[1], &z__[z_offset], &zs[1], &x[1], &diag[1], &g[1]);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L110;
	}
      g[i__] += b[i__];
    L110:
      ;
    }
  /* 
   *    ---------- 
   *    iteration 1 
   *    ------test de convergence 
   */
  s0 = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L120;
	}
      s0 += g[i__] * g[i__] / diag[i__];
    L120:
      ;
    }
  if (s0 < 1e-18)
    {
      return 0;
    }
  s1 = s0;
  /*    ------recherche de la direction de descente 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L130;
	}
      d__[i__] = -g[i__] / diag[i__];
    L130:
      ;
    }
  /* 
   *    ------step length 
   */
  dg = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L135;
	}
      dg += d__[i__] * g[i__];
    L135:
      ;
    }
  optim_calbx (n, &index[1], &indic[1], nt, np, &y[y_offset], &s[s_offset],
	       &ys[1], &z__[z_offset], &zs[1], &d__[1], &diag[1], &g[1]);
  d2a = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L140;
	}
      d2a += d__[i__] * g[i__];
    L140:
      ;
    }
  /* 
   */
  ro = -dg / d2a;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L150;
	}
      x[i__] += ro * d__[i__];
    L150:
      ;
    }
  optim_calbx (n, &index[1], &indic[1], nt, np, &y[y_offset], &s[s_offset],
	       &ys[1], &z__[z_offset], &zs[1], &x[1], &diag[1], &g[1]);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L170;
	}
      g[i__] += b[i__];
    L170:
      ;
    }
  /* 
   *    iteration k 
   */
  iter = 0;
  itmax = *np << 1;
 L10:
  ++iter;
  if (iter > itmax)
    {
      return 0;
    }
  /*    ------test de convergence 
   */
  s2 = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L200;
	}
      s2 += g[i__] * g[i__] / diag[i__];
    L200:
      ;
    }
  if (s2 / s0 < *eps)
    {
      return 0;
    }
  /*    ------recherche de la direction de descente 
   */
  beta = s2 / s1;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L210;
	}
      d__[i__] = -g[i__] / diag[i__] + beta * d__[i__];
    L210:
      ;
    }
  s1 = s2;
  /* 
   *    -----step length 
   */
  dg = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L215;
	}
      dg += d__[i__] * g[i__];
    L215:
      ;
    }
  optim_calbx (n, &index[1], &indic[1], nt, np, &y[y_offset], &s[s_offset],
	       &ys[1], &z__[z_offset], &zs[1], &d__[1], &diag[1], &g[1]);
  d2a = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L220;
	}
      d2a += d__[i__] * g[i__];
    L220:
      ;
    }
  /* 
   */
  ro = -dg / d2a;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L230;
	}
      x[i__] += ro * d__[i__];
    L230:
      ;
    }
  optim_calbx (n, &index[1], &indic[1], nt, np, &y[y_offset], &s[s_offset],
	       &ys[1], &z__[z_offset], &zs[1], &x[1], &diag[1], &g[1]);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L240;
	}
      g[i__] += b[i__];
    L240:
      ;
    }
  goto L10;
}				/* gcp_ */


/* 
 *    fonction : {bx}=[b]*{x}. 
 *               [b] est definie par les vecteurs 
 *               ({y}(i),{s}(i),{z}(i), i=1,np) et {diag} 
 * 
 *    Copyright INRIA 
 * 
 * 
 */

static int optim_calbx (int *n, int *index, int *indic, int *nt, int *np, double *y,
			double *s, double *ys, double *z__, double *zs, double *x,
			double *diag, double *bx)
{
  /* System generated locals */
  int y_dim1, y_offset, s_dim1, s_offset, z_dim1, z_offset, i__1, i__2;

  /* Local variables */
  int i__, j, ii;
  double yx, zx;

  /* Parameter adjustments */
  --bx;
  --diag;
  --x;
  --indic;
  --zs;
  z_dim1 = *nt;
  z_offset = z_dim1 + 1;
  z__ -= z_offset;
  --ys;
  s_dim1 = *nt;
  s_offset = s_dim1 + 1;
  s -= s_offset;
  y_dim1 = *nt;
  y_offset = y_dim1 + 1;
  y -= y_offset;
  --index;

  /* Function Body */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > 0)
	{
	  goto L100;
	}
      bx[i__] = diag[i__] * x[i__];
    L100:
      ;
    }
  /* 
   */
  i__1 = *np;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ii = index[i__];
      yx = 0.;
      zx = 0.;
      i__2 = *n;
      for (j = 1; j <= i__2; ++j)
	{
	  if (indic[j] > 0)
	    {
	      goto L120;
	    }
	  yx += y[ii + j * y_dim1] * x[j];
	  zx += z__[ii + j * z_dim1] * x[j];
	L120:
	  ;
	}
      /* 
       */
      i__2 = *n;
      for (j = 1; j <= i__2; ++j)
	{
	  if (indic[j] > 0)
	    {
	      goto L130;
	    }
	  bx[j] =
	    bx[j] + yx * y[ii + j * y_dim1] / ys[ii] - zx * z__[ii +
								j * z_dim1] /
	    zs[ii];
	L130:
	  ;
	}
    }
  return 0;
}



static double optim_rednor (int *n, double *binf, double *bsup, double *x,
			    double *epsx, double *g)
{
  int i__1, i__;
  double ret_val, aa;

  --g;
  --epsx;
  --x;
  --bsup;
  --binf;

  ret_val = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      aa = g[i__];
      if (x[i__] - binf[i__] <= epsx[i__])
	{
	  aa = Min (0., aa);
	}
      if (bsup[i__] - x[i__] <= epsx[i__])
	{
	  aa = Max (0., aa);
	}
      ret_val += aa * aa;
    }
  ret_val = sqrt (ret_val);
  return ret_val;
}	



static int optim_relvar (int *ind, int *n, double *x, double *binf, double *bsup,
			 double *x2, double *g, double *diag, int *imp, int *io,
			 int *ibloc, int *izag, int *iter, int *nfac, int *irit)
{
  int i__1;
  double d__1;

  /* Local variables */
  int ifac;
  double frac;
  int izag1, idfac, i__, k;
  double d1, d2, dd, bi, bs, ep, eps1;

  /*    Copyright INRIA 
   *    determination des variables a relacher par meth bertsekas 
   *    x2 vect de travail de dim n 
   *    ind: =1 si relachement des vars 
   *         =0 sinon 
   * 
   *    calcul eps1 
   */
  /* Parameter adjustments */
  --ibloc;
  --diag;
  --g;
  --x2;
  --bsup;
  --binf;
  --x;

  /* Function Body */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L10: */
      x2[i__] = x[i__] - (d__1 = g[i__], Abs (d__1)) * g[i__] / diag[i__];
    }
  optim_proj (n, &binf[1], &bsup[1], &x2[1]);
  eps1 = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L20: */
      eps1 += (d__1 = x2[i__] - x[i__], Abs (d__1));
    }
  if (*imp > 2)
    {
      Sciprintf( "relvar1. valeur de eps1=%15.7f\n", eps1);
    }
  /*    nfac nombre de lignes factorisees (nr pour ajournd) 
   */
  ifac = 0;
  idfac = 0;
  k = 0;
  frac = .10000000000000001;
  i__1 = *n;
  for (k = 1; k <= i__1; ++k)
    {
      bi = binf[k];
      bs = bsup[k];
      d1 = x[k] - bi;
      d2 = bs - x[k];
      dd = (bs - bi) * frac;
      ep = Min (eps1, dd);
      if (d1 > ep)
	{
	  goto L324;
	}
      if (g[k] > 0.)
	{
	  goto L330;
	}
      goto L335;
    L324:
      if (d2 > ep)
	{
	  goto L335;
	}
      if (g[k] > 0.)
	{
	  goto L335;
	}
      goto L330;
      /*    on defactorise si necessaire 
       */
    L330:
      if (ibloc[k] > 0)
	{
	  goto L340;
	}
      ibloc[k] = *iter;
      ++idfac;
      --(*nfac);
      *ind = 1;
      if (*imp >= 4)
	{
	  Sciprintf("relvar1: defactorisation de x(%d)=%15.7f\n",    k,    x[k]);
	}
      goto L340;
      /*    on factorise 
       */
    L335:
      if (*irit == 0)
	{
	  goto L340;
	}
      if (ibloc[k] <= 0)
	{
	  goto L340;
	}
      izag1 = *iter - ibloc[k];
      if (*izag >= izag1)
	{
	  goto L340;
	}
      ++ifac;
      ++(*nfac);
      ibloc[k] = -(*iter);
      if (*imp >= 4)
	{
	  Sciprintf( "relvar1: on factorise l indice %d\n",    k);
	}
    L340:
      ;
    }
  if (*imp >= 2 && (ifac > 0 || idfac > 0))
    {
      Sciprintf("relvar1. nbre fact %d, nbre defact %d nbre var factorisees %d\n",
		ifac,
		idfac,
		(*nfac));
    }
  *ind = 1;
  if (ifac == 0 && idfac == 0)
    {
      *ind = 0;
    }
  return 0;
}

/*    mise a l echelle de diag par la methode de shanno-phua 
 *    calcul du facteur d echelle scal 
 *     diag=(y,(diag-1)y)/(y,s)*diag 
 * 
 */


static int optim_shanph (double *diag, int *n, int *nt, int *np, double *y, double *s,
			 double *ys, double *scal, int *index, int *io, int *imp)
{
  int y_dim1, y_offset, s_dim1, s_offset, i__1;
  double d__1;
  int i__;
  double cof;
  int inp;

  /* Parameter adjustments */
  --diag;
  --index;
  --ys;
  s_dim1 = *nt;
  s_offset = s_dim1 + 1;
  s -= s_offset;
  y_dim1 = *nt;
  y_offset = y_dim1 + 1;
  y -= y_offset;

  /* Function Body */
  inp = index[*np];
  cof = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L203: */
      /*Computing 2nd power 
       */
      d__1 = y[inp + i__ * y_dim1];
      cof += d__1 * d__1 / diag[i__];
    }
  cof /= ys[inp];
  if (*imp > 3)
    {
      Sciprintf("gcbd. facteur d echelle=%15.7f\n",
		cof);
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L205: */
      diag[i__] = cof * diag[i__];
    }
  *scal = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L206: */
      *scal += diag[i__];
    }
  *scal = *n / *scal;
  return 0;
}				/* shanph_ */




/*    mise a jour de diag par la methode de bfgs diagonal 
 *    utiliser a la suite de la correction de powell 
 *    condm borne sup du conditionnement de diag 
 *    param borne inf rapport reduction diag(i) 
 * 
 *    Copyright INRIA 
 * 
 * 
 */

static int optim_bfgsd (double *diag, int *n, int *nt, int *np, double *y, double *s,
			double *ys, double *condm, double *param, double *zero,
			int *index)
{
  /* System generated locals */
  int y_dim1, y_offset, s_dim1, s_offset, i__1;
  double d__1, d__2;

  /* Local variables */
  double dmin__, omeg, dmax__;
  int i__;
  double dd, dd1, ys1;
  int inp;
  double sds, sds1;

  /* Parameter adjustments */
  --diag;
  --index;
  --ys;
  s_dim1 = *nt;
  s_offset = s_dim1 + 1;
  s -= s_offset;
  y_dim1 = *nt;
  y_offset = y_dim1 + 1;
  y -= y_offset;

  /* Function Body */
  inp = index[*np];
  ys1 = 1. / ys[inp];
  sds = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L10: */
      /*Computing 2nd power 
       */
      d__1 = s[inp + i__ * s_dim1];
      sds += diag[i__] * (d__1 * d__1);
    }
  sds1 = 1. / sds;
  dmin__ = 1e25;
  dmax__ = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      dd1 = *param * diag[i__];
      dd1 += *zero * 1e3;
      /*Computing 2nd power 
       */
      d__1 = y[inp + i__ * y_dim1];
      /*Computing 2nd power 
       */
      d__2 = diag[i__] * s[inp + i__ * s_dim1];
      dd = diag[i__] + ys1 * (d__1 * d__1) - sds1 * (d__2 * d__2);
      diag[i__] = dd;
      /*    sauvegarde positivite 
       */
      if (dd <= dd1)
	{
	  diag[i__] = dd1;
	}
      /*    calcul conditionnement 
       */
      if (diag[i__] < dmin__)
	{
	  dmin__ = diag[i__];
	}
      if (diag[i__] > dmax__)
	{
	  dmax__ = diag[i__];
	}
      /* L20: */
    }
  /*    limitation du conditionnement 
   */
  if (*condm * dmin__ / dmax__ > 1.)
    {
      return 0;
    }
  omeg = log (*condm) / log (dmax__ / dmin__);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L30: */
      diag[i__] = pow(diag[i__],omeg);
    }
  return 0;
} 



/* 
 *    mise a jour de ({z}(i),zs(i), i=1,np). 
 *    {z}(i)=[b](i-1)*{s}(i), [b](i) est definie par ({y}(j),{s}(j),{z}(j) 
 *    , j=1,i) et {diag}. 
 *    zs(i)=<z>(i)*{s}(i) 
 * 
 *    Copyright INRIA 
 * 
 * 
 */

static int  optim_majz (int *n, int *np, int *nt, double *y, double *s, double *z__,
			double *ys, double *zs, double *diag, int *index)
{
  /* System generated locals */
  int y_dim1, y_offset, s_dim1, s_offset, z_dim1, z_offset, i__1, i__2, i__3;

  /* Local variables */
  int i__, j, l, jj, jl;
  double psy, psz;

  /* Parameter adjustments */
  --diag;
  --index;
  --zs;
  --ys;
  z_dim1 = *nt;
  z_offset = z_dim1 + 1;
  z__ -= z_offset;
  s_dim1 = *nt;
  s_offset = s_dim1 + 1;
  s -= s_offset;
  y_dim1 = *nt;
  y_offset = y_dim1 + 1;
  y -= y_offset;

  /* Function Body */
  l = index[1];
  i__1 = *n;
  for (jj = 1; jj <= i__1; ++jj)
    {
      z__[l + jj * z_dim1] = diag[jj] * s[l + jj * s_dim1];
      /* L100: */
    }
  /* 
   */
  zs[l] = 0.;
  i__1 = *n;
  for (jj = 1; jj <= i__1; ++jj)
    {
      zs[l] += z__[l + jj * z_dim1] * s[l + jj * s_dim1];
      /* L110: */
    }
  /* 
   * 
   */
  if (*np == 1)
    {
      return 0;
    }
  /* 
   */
  i__1 = *np;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      l = index[i__];
      i__2 = *n;
      for (jj = 1; jj <= i__2; ++jj)
	{
	  z__[l + jj * z_dim1] = diag[jj] * s[l + jj * s_dim1];
	  /* L210: */
	}
      i__2 = i__ - 1;
      for (j = 1; j <= i__2; ++j)
	{
	  psy = 0.;
	  psz = 0.;
	  jl = index[j];
	  i__3 = *n;
	  for (jj = 1; jj <= i__3; ++jj)
	    {
	      psy += y[jl + jj * y_dim1] * s[l + jj * s_dim1];
	      psz += z__[jl + jj * z_dim1] * s[l + jj * s_dim1];
	      /* L230: */
	    }
	  i__3 = *n;
	  for (jj = 1; jj <= i__3; ++jj)
	    {
	      z__[l + jj * z_dim1] =
		z__[l + jj * z_dim1] + psy * y[jl + jj * y_dim1] / ys[jl] -
		psz * z__[jl + jj * z_dim1] / zs[jl];
	      /* L240: */
	    }
	  /* L220: */
	}
      /* 
       */
      zs[l] = 0.;
      i__2 = *n;
      for (jj = 1; jj <= i__2; ++jj)
	{
	  zs[l] += z__[l + jj * z_dim1] * s[l + jj * s_dim1];
	  /* L250: */
	}
      /* L200: */
    }
  /* 
   */
  return 0;
}

