#include "optim.h"

static int optim_strang (opt_prosca prosca, int *n, int *m, double *depl, int *jmin, int *jmax,
			 double *precon, double *alpha, double *ybar, double *sbar, opt_simul_data *optim_data);

static int optim_n1qn2a (opt_simul simul, opt_prosca prosca, int *n, double *x, double *f,
			 double *g, double *dxmin, double *df1, double *epsg,
			 int *impres, int *io, int *mode, int *niter, int *nsim, int *m,
			 double *d__, double *gg, double *aux, double *alpha,
			 double *ybar, double *sbar, opt_simul_data *optim_data);

/*!But 
 *    Minimisation sans contrainte par un algorithme 
 *    de quasi-Newton a memoire limitee. 
 *!Origine 
 *    Version 1.0 de n1qn2 (Modulopt, INRIA), septembre 1988. 
 *    Jean Charles Gilbert, Claude Lemarechal. 
 *    Copyright INRIA 
 *!Commentaires 
 *    Ce code est en principe destine aux problemes de grande taille, 
 *    n grand, mais convient egalement pour n quelconque. La methode 
 *    utilisee est du type quasi-Newton (BFGS) a encombrement variable, 
 *    ce qui permet d'utiliser au maximum la memoire declaree dispo- 
 *    nible. On estine que plus la memoire utilisee est importante, plus 
 *    rapide sera la decroissance du critere f. 
 *!Sous-routines appelees 
 *    n1qn2:    routine-chapeau qui structure la memoire declaree dispo- 
 *              nible et appelle n1qn2a, 
 *    n1qn2a:   optimiseur proprement dit, 
 *    strang:   routine de calcul de la direction de descente, 
 *    nlis0:    routine de recherche lineaire. 
 * 
 *    De son cote, l'utilisateur doit fournir: 
 *    1) une routine qui appelle le module d'optimisation n1qn2, 
 *    2) une routine de simulation, appelee simul par n1qn2, qui 
 *       calcule la valeur de f et de son gradient en un point donne, 
 *    3) une routine, appelee prosca par n1qn2, qui realise le produit 
 *       scalaire de deux vecteurs, ce produit scalaire doit etre 
 *       celui utilise pour calculer le gradient de f dans simul. 
 *!Liste d'appel 
 *    subroutine n1qn2 (simul,prosca,n,x,f,g,dxmin,df1,epsg,impres,io, 
 *   /                  mode,niter,nsim,dz,ndz,izs,rzs,dzs) 
 * 
 *    Dans la description des arguments qui suit, (e) signifie que 
 *    l'argument doit etre initialise avant l'appel de n1qn2, (s) 
 *    signifie que l'argument est une variable n'ayant de signification 
 *    qu'en sortie et (es) = (e)+(s). Les arguments du type (s) et (es) 
 *    sont en general modifies par n1qn2 et ne peuvent donc pas etre 
 *    des constantes. 
 * 
 *    simul:     Nom d'appel de la sous-routine de simulation qui 
 *               qui calcule la valeur de f et de son gradient g 
 *               a l'itere courant. Ce module doit se presenter comme 
 *               suit: 
 *                  subroutine simul (indic,n,x,f,g,izs,rzs,dzs). 
 *               Le nom de la sous-routine doit etre declare external 
 *               dans le module appelant n1qn2. Les arguments n, x, f, 
 *               g, izs, rzs et dzs ont la meme signification que ci- 
 *               dessous. N1qn2 appelle simul soit avec indic=1, dans 
 *               ce cas le simulateur fera ce qu'il veut mais ne chan- 
 *               gera pas la valeur des arguments, ou avec indic=4, dans 
 *               cas le simulateur calculera a la fois f et g. 
 *    prosca:    Nom d'appel de la sous-routine effectuant le produit 
 *               scalaire de deux vecteurs u et v. Ce module doit se 
 *               presente sous la forme: 
 *                  subroutine prosca (n,u,v,ps,izs,rzs,dzs). 
 *               Le nom de la sous-routine doit etre declare external 
 *               dans le module appelant n1qn2. Les argument n, izs, 
 *               rzs et dzs ont la meme signification que ci-dessous. 
 *               Les arguments u, v et ps sont des vecteurs de dimension 
 *               n du type double precision. Ps donne le produit 
 *               scalaire de u et v. 
 *    n(e):      Scalaire du type int. Donne la dimension n de la 
 *               variable x. 
 *    x(es):     Vecteur de dimension n du type double precision. En 
 *               entree il faut fournir la valeur du point initial, en 
 *               sortie, c'est le point final calcule par n1qn2. 
 *    f(es):     Scalaire du type double precision. En entree, c'est la 
 *               valeur de f en x (initial), valeur que l'on obtiendra 
 *               en appelant le simulateur simul avant d'appeler n1qn2. 
 *               En sortie, c'est la valeur de f en x (final). 
 *    g(es):     Vecteur de dimension n du type double precision. 
 *               En entree, il faut fournir la valeur du gradient de f 
 *               en x (initial), valeur que l'on obtiendra en appelant 
 *               le simulateur simul avant d'appeler n1qn2. En sortie en 
 *               mode 1, c'est la valeur du gradient de f en x (final). 
 *    dxmin(e):  Scalaire du type double precision, strictement positif. 
 *               Cet argument definit la resolution sur x en norme 
 *               l-infini: deux points dont la distance en norme l- 
 *               infini est superieure a dxmin seront consideres comme 
 *               non distinguables par la routine de recherche lineaire. 
 *    df1(e):    Scalaire du type double precision, strictement positif. 
 *               Cet argument donne une estimation de la diminution 
 *               escomptee pour f lors de la premiere iteration. 
 *    epsg(es):  Scalaire du type double precision, strictement positif 
 *               et strictement inferieur a 1. Sa valeur en entree, 
 *               determine le test d'arret qui porte sur la norme 
 *               (prosca) du gradient. Le minimiseur considere que la 
 *               convergence est atteinte en x(k) et s'arrete en mode 1 
 *               si E(k) := |g(k)|/|g(1)| < epsg, ou g(1) et g(k) sont 
 *               les gradients au point d'entree et a l'iteration k, 
 *               respectivement. En sortie, epsg = E(k). 
 *    impres(e): Scalaire du type int qui controle les sorties. 
 *               <0:  Rien n'est imprime et n1qn2 appelle le simulateur 
 *                    avec indic=1 toutes les (-impres) iterations. 
 *               =0:  Rien n'est imprime. 
 *               >=1: Impressions initiales et finales, messages 
 *                    d'erreurs. 
 *               >=3: Une ligne d'impression par iteration donnant 
 *                    l'ordre k de l'iteration courante menant de x(k) 
 *                    a x(k+1), le nombre d'appels au simulateur avant 
 *                    cette iteration, la valeur du critere f et sa 
 *                    derivee directionnelle suivant d(k). 
 *               >=4: Impression de nlis0. 
 *               >=5: Impressions supplentaires en fin d'iteration k: 
 *                    le test d'arret E(k+1), prosca(y(k),s(k)) qui 
 *                    doit etre positif, le facteur de Oren-Spedicato 
 *                    et l'angle de la direction de descente d(k) avec 
 *                    -g(k). 
 *    io(e):     Scalaire du type int qui sera pris comme numero 
 *               de canal de sortie pour les impressions controlees 
 *               par impres. 
 *    mode(s):   Scalaire du type int donnant le mode de sortie de 
 *               n1qn2. 
 *               <0: Impossibilite de poursuivre la recherche lineaire 
 *                   car le simulateur a repondu avec indic<0. Mode 
 *                   renvoie cette valeur de indic. 
 *               =0: Arret demande par le simulateur qui a repondu avec 
 *                   indic=0. 
 *               =1: Fin normale avec test sur le gradient satisfait. 
 *               =2: Arguments d'entree mal initialises. Il peut s'agir 
 *                   de n<=0, niter<=0, nsim<=0, dxmin<=0, epsg<=0, 
 *                   epsg>1 ou de nrz<5n+1 (pas assez de memoire 
 *                   allouee). 
 *               =3: La recherche lineaire a ete bloquee sur tmax=10^20 
 *                   (mode tres peu probable). 
 *               =4: Nombre maximal d'iterations atteint. 
 *               =5: Nombre maximal de simulations atteint. 
 *               =6: Arret sur dxmin lors de la recherche lineaire. Ce 
 *                   mode de sortie peut avoir des origines tres 
 *                   diverses. Si le nombre d'essais de pas lors de la 
 *                   derniere recherche lineaire est faible, cela peut 
 *                   signifier que dxmin a ete pris trop grand. Il peut 
 *                   aussi s'agir d'erreurs ou d'imprecision dans le 
 *                   calcul du gradient. Dans ce cas, la direction de 
 *                   recherche d(k) peut ne plus etre une direction de 
 *                   descente de f en x(k), etant donne que n1qn2 
 *                   s'autorise des directions d(k) pouvant faire avec 
 *                   -g(k) un angle proche de 90 degres. On veillera 
 *                   donc a calculer le gradient avec precision. 
 *               =7: Soit (g,d) soit (y,s) ne sont pas positifs (mode 
 *                   de sortie tres improbable). 
 *    niter(es): Scalaire du type int, strictement positif. En 
 *               entree, c'est le nombre maximal d'iterations admis. 
 *               En sortie, c'est le nombre d'iterations effectuees. 
 *    nsim(es):  Scalaire du type int, strictement positif. En 
 *               entree, c'est le nombre maximal de simulations admis. 
 *               En sortie, c'est le nombre de simulations effectuees. 
 *    rz(s):     Vecteur de dimension nrz du type double precision. 
 *               C'est l'adresse d'une zone de travail pour n1qn2. 
 *    nrz(e):    Scalaire du type int, strictement positif, donnant 
 *               la dimension de la zone de travail rz. En plus des 
 *               vecteurs x et g donnes en arguments, n1qn2 a besoin 
 *               d'une zone de travail composee d'au moins trois 
 *               vecteurs de dimension n et chaque mise a jour demandee 
 *               necessite 1 scalaire et 2 vecteurs de dimension n 
 *               supplementaires. Donc si m est le nombre de mises a 
 *               jour desire pour la construction de la metrique 
 *               locale, il faudra prendre 
 *                  nrz >= 3*n + m*(2*n+1). 
 *               En fait, le nombre m est determine par n1qn2 qui prend: 
 *                  m = partie entiere par defaut de ((nrz-3*n)/(2*n+1)) 
 *               Ce nombre doit etre >= 1. Il faut donc nrz >= 5*n +1, 
 *               sinon n1qn2 s'arrete en mode 2. 
 *    izs, rzs, dzs:  Adresses de zones-memoire respectivement du type 
 *               int, real et double precision. Elles sont reservees 
 *               a l'utilisateur. N1qn2 ne les utilise pas et les 
 *               transmet a simul et prosca. 
 */




int optim_n1qn2 (opt_simul simul, opt_prosca prosca, int *n, double *x, double *f, double *g,
		 double *dxmin, double *df1, double *epsg, int *impres, int *io,
		 int *mode, int *niter, int *nsim, double *dz, int *ndz, opt_simul_data *optim_data)
{
  /* Local variables */
  int iaux, ndzu, m, isbar;
  int iybar;
  double r1, r2;
  int l1memo, id, ialpha;
  double ps;
  int igg;


  /* Parameter adjustments */
  --dz;
  --g;
  --x;

  /* Function Body */
  if (*impres >= 1)
    {
      
      Sciprintf("n1qn2: point d'entree\n\tdimension du probleme (n)              :,%d\n\tprecision absolue en x (dxmin)         :,%9.2f\n\tdecroissance attendue pour f (df1)     :,%9.2f\n\tprecision relative en g (epsg)         :,%9.2f\n\tnombre maximal d'iterations (niter)    :,%d\n\tnombre maximal d'appels a simul (nsim) :,%d\n\tniveau d'impression (impres)           :,%d\n",
		(*n),
		(*dxmin),
		(*df1),
		(*epsg),
		(*niter),
		(*nsim),
		(*impres));
    }
  if (*n <= 0 || *niter <= 0 || *nsim <= 0 || *dxmin <= 0. || *epsg <= 0.
      || *epsg > 1.)
    {
      *mode = 2;
      if (*impres >= 1)
	{
	  Sciprintf("n1qn2 : appel incoherent\n");
	}
      goto L904;
    }
  if (*ndz < *n * 5 + 1)
    {
      *mode = 2;
      if (*impres >= 1)
	{
	  Sciprintf("n1qn2: memoire allouee insuffisante\n");
	}
      goto L904;
    }
  /* 
   *---- calcul de m et des pointeurs subdivisant la zone de travail dz 
   * 
   */
  ndzu = *ndz - *n * 3;
  l1memo = (*n << 1) + 1;
  m = ndzu / l1memo;
  ndzu = m * l1memo + *n * 3;
  if (*impres >= 1)
    {
      Sciprintf("n1qn2: memoire allouee (ndz) %d, memoire utilisee: %d, nombre de mises a jour: %d\n",
		(*ndz),
		ndzu,
		m);
    }
  id = 1;
  igg = id + *n;
  iaux = igg + *n;
  ialpha = iaux + *n;
  iybar = ialpha + m;
  isbar = iybar + *n * m;
  /* 
   *---- appel du code d'optimisation 
   * 
   */
  optim_n1qn2a ( simul,  prosca, n, &x[1], f, &g[1], dxmin, df1,
		 epsg, impres, io, mode, niter, nsim, &m, &dz[id], &dz[igg],
		 &dz[iaux], &dz[ialpha], &dz[iybar], &dz[isbar],optim_data);
  /* 
   *---- impressions finales 
   * 
   */
 L904:
  if (*impres >= 1)
    {
      Sciprintf("n1qn2: sortie en mode %d, nombre d'iterations %d nombre d'appels a simul %d precision relative atteinte sur g: %9.2f\n",
		(*mode),
		(*niter),
		(*nsim),
		(*epsg));
    }
  (*prosca) (n, &x[1], &x[1], &ps,optim_data);
  r1 = sqrt (ps);
  (*prosca) (n, &g[1], &g[1], &ps,optim_data);
  r2 = sqrt (ps);
  if (*impres >= 1)
    {
      Sciprintf("n1qn2: norme de x =%15.8f, f=%15.8f, norme de g=%15.8f\n",
		r1,
		(*f),
		r2);
    }
  /* 
   */
  return 0;
}	






static int optim_n1qn2a (opt_simul simul, opt_prosca prosca, int *n, double *x, double *f,
			 double *g, double *dxmin, double *df1, double *epsg,
			 int *impres, int *io, int *mode, int *niter, int *nsim, int *m,
			 double *d__, double *gg, double *aux, double *alpha,
			 double *ybar, double *sbar, opt_simul_data *optim_data)
{
  double c_b25 = .9;
  double c_b26 = 1e-4;

  /* System generated locals */
  int ybar_dim1, ybar_offset, sbar_dim1, sbar_offset, i__1;
  double d__1, d__2, d__3;

  /* Local variables */
  int jmin, jmax, isim, iter;
  double tmin, tmax;
  int i__;
  double r__, t;
  int indic;
  double gnorm, d1, ff, ps;
  int moderl;
  double precon;
  double hp0, ps2, eps1;

  /*    Copyright INRIA 
   *---- 
   * 
   *    code d'optimisation proprement dit 
   * 
   *---- 
   * 
   *    arguments 
   * 
   * 
   *    variables locales 
   * 
   * 
   *---- parametres 
   * 
   * 
   *---- initialisation 
   * 
   */
  /* Parameter adjustments */
  --aux;
  --gg;
  --d__;
  --g;
  --x;
  sbar_dim1 = *n;
  sbar_offset = sbar_dim1 + 1;
  sbar -= sbar_offset;
  ybar_dim1 = *n;
  ybar_offset = ybar_dim1 + 1;
  ybar -= ybar_offset;
  --alpha;

  /* Function Body */
  iter = 0;
  isim = 1;
  (*prosca) (n, &g[1], &g[1], &ps,optim_data);
  gnorm = sqrt (ps);
  if (*impres >= 1)
    {
      Sciprintf("n1qn2: f=%15.8f norme de g=%15.8f\n",
		(*f),
		gnorm);
    }
  /* 
   *    ---- direction de descente initiale 
   *         (avec mise a l'echelle) 
   * 
   *Computing 2nd power 
   */
  d__1 = gnorm;
  precon = *df1 * 2. / (d__1 * d__1);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      d__[i__] = -g[i__] * precon;
      /* L10: */
    }
  if (*impres >= 5)
    {
      Sciprintf("n1qn2: direction de descente -g: precon =%10.3f\n",precon);
    }
  if (*impres == 3)
    {
      Sciprintf("n1qn2:\n");
    }
  if (*impres == 4)
    {
      Sciprintf( "n1qn2:\n");
    }
  /* 
   *    ---- initialisation pour nlis0 
   * 
   */
  tmax = 1e20;
  (*prosca) (n, &d__[1], &g[1], &hp0,optim_data);
  /* 
   *    ---- initialisation pour strang 
   * 
   */
  jmin = 1;
  jmax = 0;
  /* 
   *---- debut de l'iteration. on cherche x(k+1) de la forme x(k) + t*d, 
   *    avec t > 0. on connait d. 
   * 
   *    debut de la boucle: etiquette 100, 
   *    sortie de la boucle: goto 1000. 
   * 
   */
 L100:
  ++iter;
  if (*impres < 0)
    {
      if (iter % (-(*impres)) == 0)
	{
	  indic = 1;
	  (*simul) (&indic, n, &x[1], f, &g[1],optim_data);
	  goto L100;
	}
    }
  if (*impres >= 5)
    {
      Sciprintf( "n1qn2:\n");
    }
  if (*impres >= 4)
    { 
      Sciprintf( "n1qn2:\n");
    }
  if (*impres >= 3)
    {
      Sciprintf("n1qn2: iter=%d, simul=%d, f=%15.8f,  h'(0)=%12.5f\n",
		iter,
		isim,
		(*f),
		hp0);
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      gg[i__] = g[i__];
      /* L101: */
    }
  ff = *f;
  /* 
   *    ---- recherche lineaire et nouveau point x(k+1) 
   * 
   */
  if (*impres >= 5)
    {
      Sciprintf("n1qn2: recherche lineaire\n");
    }
  /* 
   *        ---- calcul de tmin 
   * 
   */
  tmin = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /*Computing MAX 
       */
      d__2 = tmin, d__3 = (d__1 = d__[i__], Abs (d__1));
      tmin = Max (d__2, d__3);
      /* L200: */
    }
  tmin = *dxmin / tmin;
  t = 1.;
  d1 = hp0;
  /* 
   */
  optim_nlis0 (n,simul, prosca, &x[1], f, &d1, &t, &tmin,
	       &tmax, &d__[1], &g[1], &c_b25, &c_b26, impres, io, &moderl,
	       &isim, nsim, &aux[1], optim_data);
  /* 
   *        ---- nlis0 renvoie les nouvelles valeurs de x, f et g 
   * 
   */
  if (moderl != 0)
    {
      if (moderl < 0)
	{
	  /* 
	   *            ---- calcul impossible 
	   *                 t, g: ou les calculs sont impossible 
	   *                 x, f: ceux du t_gauche (donc f <= ff) 
	   * 
	   */
	  *mode = moderl;
	}
      else if (moderl == 1)
	{
	  /* 
	   *            ---- descente bloquee sur tmax 
	   *                 [sortie rare (!!) d'apres le code de nlis0] 
	   * 
	   */
	  *mode = 3;
	  if (*impres >= 1)
	    {
	      Sciprintf("n1qn2: (iteration %d): recherche lineaire bloquee sur tmax: reduire l'echelle\n",
			iter);
	    }
	}
      else if (moderl == 4)
	{
	  /* 
	   *            ---- nsim atteint 
	   *                 x, f: ceux du t_gauche (donc f <= ff) 
	   * 
	   */
	  *mode = 5;
	}
      else if (moderl == 5)
	{
	  /* 
	   *            ---- arret demande par l'utilisateur (indic = 0) 
	   *                 x, f: ceux en sortie du simulateur 
	   * 
	   */
	  *mode = 0;
	}
      else if (moderl == 6)
	{
	  /* 
	   *            ---- arret sur dxmin ou appel incoherent 
	   *                 x, f: ceux du t_gauche (donc f <= ff) 
	   * 
	   */
	  *mode = 6;
	}
      goto L1000;
    }
  /* 
   *    ---- tests d'arret 
   * 
   */
  (*prosca) (n, &g[1], &g[1], &ps,optim_data);
  eps1 = sqrt (ps) / gnorm;
  /* 
   */
  if (*impres >= 5)
    {
      Sciprintf( "n1qn2: test d'arret sur g: %12.5f\n", eps1);
    }
  if (eps1 < *epsg)
    {
      *mode = 1;
      goto L1000;
    }
  if (iter == *niter)
    {
      *mode = 4;
      if (*impres >= 1)
	{
	  Sciprintf("n1qn2 (iteration %d): nombre maximal d'iterations atteint\n",
		    iter);
	}
      goto L1000;
    }
  if (isim >= *nsim)
    {
      *mode = 5;
      if (*impres >= 1)
	{
	  Sciprintf("n1qn2 (iteration %d): %d appels a simul (nombre maximal atteint)\n",
		    iter,
		    isim);
	}
      goto L1000;
    }
  /* 
   *    ---- mise a jour de la matrice 
   * 
   */
  ++jmax;
  if (iter > *m)
    {
      ++jmin;
      if (jmin > *m)
	{
	  jmin -= *m;
	}
      if (jmax > *m)
	{
	  jmax -= *m;
	}
    }
  /* 
   *        ---- y, s et (y,s) 
   * 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      sbar[i__ + jmax * sbar_dim1] = t * d__[i__];
      ybar[i__ + jmax * ybar_dim1] = g[i__] - gg[i__];
      /* L400: */
    }
  (*prosca) (n, &ybar[jmax * ybar_dim1 + 1], &sbar[jmax * sbar_dim1 + 1], &d1,optim_data);
  if (d1 <= 0.)
    {
      *mode = 7;
      if (*impres >= 1)
	{
	  Sciprintf("n1qn2 (iteration %d): le produit scalaire (y,s) =%12.5f n'est pas positif\n",
		    iter,
		    d1);
	}
      goto L1000;
    }
  /* 
   *        ---- precon: facteur de mise a l'echelle 
   * 
   */
  (*prosca) (n, &ybar[jmax * ybar_dim1 + 1], &ybar[jmax * ybar_dim1 + 1], &ps,optim_data);
  precon = d1 / ps;
  if (*impres >= 5)
    {
      Sciprintf("n1qn2: mise a jour: (y,s) =%10.3f, Oren-Spedicato =%10.3f\n",
		d1,
		precon);
    }
  /* 
   *        ---- ybar, sbar 
   * 
   */
  d1 = sqrt (1. / d1);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      sbar[i__ + jmax * sbar_dim1] = d1 * sbar[i__ + jmax * sbar_dim1];
      ybar[i__ + jmax * ybar_dim1] = d1 * ybar[i__ + jmax * ybar_dim1];
      /* L410: */
    }
  /* 
   *    ---- calcul de la nouvelle direction de descente d = - h.g 
   * 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      d__[i__] = -g[i__];
      /* L510: */
    }
  optim_strang ( prosca, n, m, &d__[1], &jmin, &jmax, &precon,
		 &alpha[1], &ybar[ybar_offset], &sbar[sbar_offset],optim_data);
  /* 
   *        ---- test: la direction d est-elle de descente ? 
   *            hp0 sera utilise par nlis0 
   * 
   */
  (*prosca) (n, &d__[1], &g[1], &hp0,optim_data);
  if (hp0 >= 0.)
    {
      *mode = 7;
      if (*impres >= 1)
	{
	  Sciprintf("n1qn2 (iteration %d) la direction de recherche d n'est pas de descente: (g,d)=%12.5f\n",
		    iter,
		    hp0);
	}
      goto L1000;
    }
  if (*impres >= 5)
    {
      (*prosca) (n, &g[1], &g[1], &ps,optim_data);
      ps = sqrt (ps);
      (*prosca) (n, &d__[1], &d__[1], &ps2,optim_data);
      ps2 = sqrt (ps2);
      ps = hp0 / ps / ps2;
      /*Computing MIN 
       */
      d__1 = -ps;
      ps = Min (d__1, 1.);
      ps = acos (ps);
      r__ = ps * 180. / 3.1415927;
      Sciprintf("n1qn2: direction de descente d: angle(-g,d) =%5.1f degres\n",r__);
    }
  /* 
   *---- on poursuit les iterations 
   * 
   */
  goto L100;
  /* 
   *        retour 
   * 
   */
 L1000:
  *epsg = eps1;
  *niter = iter;
  *nsim = isim;
  return 0;
}				/* n1qn2a_ */




static int optim_strang (opt_prosca prosca, int *n, int *m, double *depl, int *jmin, int *jmax,
			 double *precon, double *alpha, double *ybar, double *sbar, opt_simul_data *optim_data)
{
  /* System generated locals */
  int ybar_dim1, ybar_offset, sbar_dim1, sbar_offset, i__1, i__2;

  /* Local variables */
  int jfin, i__, j;
  double r__;
  int jp;
  double ps;

  /*---- 
   * 
   *    Calcule le produit H g ou 
   *        . H est une matrice construite par la formule de bfgs inverse 
   *          a m memoires a partir de precon fois la matrice unite dans 
   *          un espace hilbertien dont le produit scalaire est donne par 
   *          prosca 
   *          (cf. J. Nocedal, math. of comp. 35/151 (1980) 773-782) 
   *        . g est un vecteur de dimension n (en general le gradient) 
   * 
   *    Le facteur precon apparait donc comme un preconditionneur 
   *    scalaire. 
   * 
   *    delp = g (en entree), = H g (en sortie) 
   * 
   *    La matrice H est memorisee par les vecteurs des tableaux 
   *    ybar, sbar et les pointeurs jmin, jmax. 
   * 
   *    alpha(m) est une zone de travail. 
   * 
   *    izs(1),rzs(1),dzs(1) sont des zones de travail pour prosca 
   * 
   *---- 
   * 
   *        arguments 
   * 
   * 
   *        variables locales 
   * 
   * 
   */
  /* Parameter adjustments */
  --depl;
  sbar_dim1 = *n;
  sbar_offset = sbar_dim1 + 1;
  sbar -= sbar_offset;
  ybar_dim1 = *n;
  ybar_offset = ybar_dim1 + 1;
  ybar -= ybar_offset;
  --alpha;

  /* Function Body */
  jfin = *jmax;
  if (jfin < *jmin)
    {
      jfin = *jmax + *m;
    }
  /* 
   *        phase de descente 
   * 
   */
  i__1 = *jmin;
  for (j = jfin; j >= i__1; --j)
    {
      jp = j;
      if (jp > *m)
	{
	  jp -= *m;
	}
      (*prosca) (n, &depl[1], &sbar[jp * sbar_dim1 + 1], &ps,optim_data);
      alpha[jp] = ps;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
	{
	  depl[i__] -= ps * ybar[i__ + jp * ybar_dim1];
	  /* L20: */
	}
      /* L100: */
    }
  /* 
   *        preconditionnement 
   * 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      depl[i__] *= *precon;
      /* L150: */
    }
  /* 
   *        remontee 
   * 
   */
  i__1 = jfin;
  for (j = *jmin; j <= i__1; ++j)
    {
      jp = j;
      if (jp > *m)
	{
	  jp -= *m;
	}
      (*prosca) (n, &depl[1], &ybar[jp * ybar_dim1 + 1], &ps,optim_data);
      r__ = alpha[jp] - ps;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
	{
	  depl[i__] += r__ * sbar[i__ + jp * sbar_dim1];
	  /* L120: */
	}
      /* L200: */
    }
  return 0;
}
