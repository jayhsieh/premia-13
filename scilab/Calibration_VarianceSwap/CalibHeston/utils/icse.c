
/*!but 
 *    Le logiciel ICSE est un outil de resolution de problemes de 
 *    CONTROLE OPTIMAL de systemes decrits par des equations 
 *    differentielles ou algebrico-differentielles, NON LINEAIRES. 
 *    Dans la mesure ou dans la methode d'integration qu'il utilise 
 *    est inconditionnellement stable,il peut aussi etre utilise pour 
 *    resoudre des problemes de controle d'EQUATIONS AUX DERIVEES 
 *    PARTIELLES DYNAMIQUES (reaction-diffusion, par exemple), ou 
 *    plus generalement pour controler des systemes raides. 
 *    Le controle se decompose en une partie independante du temps et 
 *    une partie dependante du temps.Cette structure permet a 
 *    l'utilisateur de resoudre facilement les problemes 
 *    d'IDENTIFICATION DE PARAMETRES d'un systeme dynamique par des 
 *    methodes de MOINDRES CARRES. Diverses facilites sont d'ailleurs 
 *    prevues pour ce cas. 
 *    Dans le cas d'un SYSTEME LINEAIRE,l'integration de l'etat et de 
 *    l'etat adjoint sont effectuees de maniere a tirer parti de la 
 *    linearite de l'equation. 
 * 
 *    Resolution de problemes de controle ou d'identification de 
 *     parametres de systemes dynamiques du type: 
 *         0=fi(t,y,u),i<=nea et dyj/dt=fj(t,y,u),nea<j<=ny pour t>=t0 
 *         et y(t0)=y0 , en notant ny la dimension de l'etat y et 
 *         fk la keme composante de la fonction f 
 *     On effectue un certain nombre (nex) d'experiences identiques 
 *     au cours desquelles certaines donnees sont mesurees. 
 *     Le critere a minimiser est de la forme c(tob,ytob,ob) avec: 
 *       tob(ntob)       :instants de mesure 
 *       ytob(ny,ntob)   :valeurs de l'etat aux instants de mesure 
 *       ob(nex,ntob,nob):mesures 
 *     L'equation d'etat est discretisee par la methode 
 *     de Crank-Nicolson. 
 *     Le gradient du cout est calcule en utilisant l'etat adjoint du 
 *     systeme discretise. 
 * 
 *!origine 
 *     F.Bonnans,G.Launay, INRIA, 1987 
 * 
 *! DESCRIPTION FORMELLE DES PROBLEMES DE CONTROLE CONSIDERES 
 *      L'equation du systeme est de la forme 
 *             0     =  fi(t,y(t),uc,uv(t)) ,  i<=nea, 
 *         dyi(t)/dt =  fi(t,y(t),uc,uv(t)) ,  i> nea, 
 *      et 
 *         t0<= t <=tf ,    y(t0)=y0 . 
 *    avec : 
 *      y(t )  : etat du systeme, 
 *      uc     : partie  du controle independante du temps 
 *              (controle constant), 
 *      uv     : partie  du controle dependante du temps 
 *               (controle variable) 
 *      t0, tf : instant initial et instant final, 
 *      y0     : etat initial. Il est soit fixe,soit fonction 
 *               du controle [y0=ei(uc,uv)] 
 * 
 *      Le critere a minimiser est de la forme  : 
 * 
 *                 tp 
 *                  ____ 
 *                 \ 
 *                  |    c2(ti,y(ti),uc) . 
 *                 /___ 
 *                 ti =t1 
 * 
 *       Sont resolus, soit des problemes sans contraintes,  soit des 
 *    problemes comportant des contraintes de borne sur le controle.On 
 *    peut aussi utiliser ICSE pour resoudre des problemes comportant 
 *    des contraintes sur  l'etat, en traitant  ces contraintes par 
 *    penalisation ou lagrangien augmente [Ber,82]. 
 *       Les fonctions  f,c1,c2,ei et leurs derivees  partielles sont 
 *    fournies par l'utilisateur sous forme de subroutines Fortran. 
 * 
 * 
 * 
 *! OUTILS POUR L'IDENTIFICATION DE PARAMETRES 
 *    Le probleme de l'identification de parametres (ou d'ajustement 
 *    de parametres a des mesures) a les caracteristiques suivantes : 
 *    le critere est fonction seulement de mesures faites a certains 
 *    instants et des valeurs de l'etat a ces instants.Seule la seconde 
 *    partie du cout intervient donc.Les mesures peuvent avoir ete 
 *    obtenues lors de plusieurs experiences.En general,le critere 
 *    est du type MOINDRES CARRES associe a une OBSERVATION LINEAIRE 
 *    .Autrement dit, il est de la forme 
 * 
 *               nex      ntob 
 *               ____     ____                              2 
 *          1    \        \      ||                       || 
 *          -     |        |     || obs*y(ti) - z(iex,ti) || , 
 *          2    /___     /___   ||                       || 
 *               iex=1    ti =t1 
 * 
 *    ou obs est  la  matrice  d'observation et  z(iex,ti) represente 
 *    l'ensemble  des  mesures  faites  lors  de l'experience iex,  a 
 *    l'instant ti.Dans ce cas,l'utilisateur n'a aucune modification a 
 *    a apporter  a la subroutine de calcul  de cout  de l'exemple de 
 *    demonstration. 
 *!liste d'appel: 
 *     icse(ind,nu,u,co,g,itv,rtv,dtv) 
 *     en entree: 
 * 
 *     ind        entier egal a 2,3,ou4 
 * 
 *     nu         entier 
 *                dimension du vecteur des parametres 
 * 
 *     u          double precision (nu) 
 *                vecteur des parametres 
 * 
 *     en sortie: 
 * 
 *     co         double precision 
 *                cout 
 * 
 *     g          double precision (nu) 
 *                gradient de la fonction de cout c 
 * 
 *     itv        entier (nitv) 
 *                tableau de travail entier 
 * 
 *     rtv        reel (nrtv) 
 *                tableau de travail reel 
 * 
 *     dtv        double precision (ndtv) 
 *                tableau de travail double precision 
 * 
 *!subroutines utilisees 
 *     Linpack    :dadd,daxpy,dcopy,dmmul,dnrm2,dscal,dset, 
 *                 dgefa,dgesl 
 *                :icse0,icse1,icse2,icscof,icsef,icsei 
 *     common/nird/nitv,nrtv,ndtv 
 *!utilisation 
 * 
 *    Le probleme a traiter doit etre defini par 3 routines 
 *    fortran ecrites par l'utilisateur : 
 * 
 *     -Second membre de l'equation d'etat : 
 *      icsef(indf,t,y,uc,uv,f,fy,fu,b,itu,dtu, 
 *    & t0,tf,dti,dtf,ermx,iu,nuc,nuv,ilin,nti,ntf,ny,nea, 
 *    & itmx,nex,nob,ntob,ntobi,nitu,ndtu) 
 *      Parametres d'entree : 
 *       indf     : vaut 1,2,3 suivant qu'on veut calculer f,fy,fu 
 *       t        : instant courant 
 *       y(ny)    : etat a un instant donne 
 *       uc(nuc)  : controle independant du temps 
 *       uv(nuv)  : controle dependant du temps, a l'instant t 
 *       b(ny)    : terme constant dans le cas lineaire quadratique 
 *      Parametres de sortie : 
 *        indf    : >0 si  le calcul s'est  correctement effectue,<=0 
 *                  sinon 
 *       f(ny)    : second membre 
 *       fy(ny,ny): jacobien de f par rapport a y 
 *       fu(ny,nuc+nuv) : derivee de f par rapport au controle 
 *      Tableaux de travail reserves a l'utilisateur : 
 *       itu(nitu): tableau entier 
 *       dtu(ndtu): tableau double precision 
 *      (nitu et ndtu sont initialises par le common icsez). 
 * 
 *     -Cout ponctuel : 
 *      icsec2(indc,nu,tob,obs,cof,ytob,ob,u,c2,c2y,g,yob,d,itu,dtu, 
 *    & t0,tf,dti,dtf,ermx,iu,nuc,nuv,ilin,nti,ntf,ny,nea, 
 *    & itmx,nex,nob,ntob,ntobi,nitu,ndtu) 
 *      Parametres d'entree : 
 *       indc     : 1 si on desire calculer c2,2 si on desire 
 *                  calculer c2y,c2u 
 *       tob      : instants de mesure 
 *       obs      : matrice d'observation 
 *       cof      : coefficients de ponderation du cout 
 *       ytob     : valeur de l'etat aux instants d'observation 
 *       ob       : mesures 
 *       u(nu)    : controle.Le controle variable est stocke a la 
 *                  suite du controle suite du constant. 
 *      Parametres de sortie : 
 *       indc     : comme pour icsec1 
 *       c2       : cout 
 *       c2y(ny,ntob) : derivee de c2 par rapport a y 
 *       g(nu)  : derivee de c2 par rapport a u 
 * 
 *     -Etat initial (s'il est variable) 
 *      icsei(indi,nui,ui,y0,y0ui,itu,dtu, 
 *    & t0,tf,dti,dtf,ermx,iu,nuc,nuv,ilin,nti,ntf,ny,nea, 
 *    & itmx,nex,nob,ntob,ntobi,nitu,ndtu) 
 *      Parametres d'entree : 
 *       indi     : 1 si on desire calculer y0, 2 si on  desire 
 *                    calculer y0ui 
 *       nui      : dimension du tableau ui defini ci-dessous, 
 *       ui       : partie du controle intervenant dans l'etat initial, 
 *                  determinee par iu;vaut uc,uv,ou [uc,uv]. 
 * 
 *      Parametres de sortie : 
 *       indc     : >0  si le calcul  s'est correctement effectue,<=0 
 *                  sinon, 
 *       y0       : etat initial, 
 *       y0ui     : derivee de l'etat initial par rapport au controle. 
 * 
 * 
 * 
 *!vue d'ensemble 
 *        Pour  utiliser la  subroutine icse, il faut  disposer d'un 
 *    optimiseur (code d'implementation d'un algorithme d'optimisation) 
 *    a  la norme  MODULOPT.Il  faut  ensuite  ecrire  le  programme 
 *    principal, constitue de quatre parties : 
 *      1. Initialisation des variables du common icsez, 
 *      2. Initialisation des tableaux itv et dtv et du common nird, 
 *      3. Appel de l'optimiseur, 
 *      4. Traitement des resultats. 
 * 
 * 
 *    1. INITIALISATION DU COMMON ICSEZ 
 *       common/icsez/t0,tf,dti,dtf,ermx,iu,nuc,nuv,ilin,nti,ntf,ny,nea, 
 *       itmx,nex,nob,ntob,ntobi,nitu,ndtu 
 * 
 *      Liste des variables a initialiser : 
 *      t0    : instant initial 
 *      tf    : instant final 
 *      dti   : premier pas de temps 
 *      dtf   : second pas de temps 
 *      ermx  : test d'arret absolu sur la valeur du  second membre 
 *              dans la resolution de l'equation d'etat 
 *      iu(5) : tableau  parametrant le  probleme : seuls iu(1:3) 
 *              sont utilises. 
 *           iu(1)=1 si l'etat initial depend du  controle constant 
 *                 0 sinon 
 *           iu(2)=1 si l'etat initial  depend du controle variable 
 *                 0 sinon 
 *           iu(3)=1 si le second membre depend du controle constant, 
 *                 0 sinon 
 * 
 *      nuc   : dimension du controle constant. 
 *      nuv   : dimension du controle variable a un instant donne. 
 *      ilin  : indicateur de linearite 
 *      nti   : nombre de pas de temps correspondant a dti (premier 
 *                pas de temps) 
 *      ntf   : nombre de pas  de temps correspondant a dtf (second 
 *              pas de temps) 
 *      ny    : dimension de l'etat a un instant donne 
 *      nea   : nombre d'equations algebriques (eventuellement nul) 
 *      itmx  : nombre  maximal  d'iterations dans la resolution 
 *              de l'equation d'etat discrete a un pas de temps 
 *              donne 
 *      nex   : nombre d'experiences effectuees 
 *      nob   : dimension du  vecteur des mesures  pour une 
 *              experience donnee en un instant donne 
 *      ntob  : nombre d'instants de mesure pour une experience donnee 
 *      ntobi : nombre d'instants de mesure correspondant a dti 
 *              (premier pas de temps) 
 *      nitu  : longueur de  itu,tableau de  travail entier reserve 
 *              a l'utilisateur 
 *      ndtu  : longueur de dtu,  tableau de travail  double 
 *              precision reserve a l'utilisateur 
 *      u(nu)       : parametres initiaux 
 *      y0(ny)      : etat initial 
 *      tob(ntob)   : instants de mesure 
 *      binf(nu)    : borne inferieures sur les parametres 
 *      bsup(nu)    : borne superieures sur les parametres 
 *      obs(nob,ny) : matrice d'observation 
 * 
 *    Bien noter que 
 *        nu = nuc + nuv*(nti+ntf+1), 
 *        nui= iu(1)*nuc+ui(2)*nuv*(nti+ntf+1) 
 *    et que les dimensions suivantes peuvent etre nulles : 
 *        nuc,nuv,ntf,nea. 
 * 
 * 
 *    2 INITIALISATION DES TABLEAUX ENTIER ET DOUBLE PRECISION. 
 *      Le tableau itv (entier) contient le tableau : 
 *      itu    dimension nitu      : reserve a l'utilisateur, 
 *    le reste du tableau etant reserve au systeme ICSE. 
 * 
 *    Le tableau dtv (reel double precision) contient les tableaux : 
 *      dtu    dimension  ndtu     : reserve a l'utilisateur, 
 *      y0                ny       : etat initial, 
 *      tob               ntob     : instants d'observation, 
 *      obs               nob,ny   : matrice d'observation, 
 *      ob            nex,ntob,nob : observations (mesures), 
 *      ech               nu       : coefficients de mise a l'echelle de 
 *                                   u, 
 *      cof              nob,ntob  : coefficients de ponderation du 
 *                                   cout, 
 * 
 *      Les dimensions nitu et ndtu sont passees par le common icsez, 
 *                     nitv et ndtv sont passees par le common nird. 
 * 
 * 
 *! 
 *    Copyright INRIA 
 * 
 * 
 * 
 *    lui et nui servent quand l'etat initial depend du controle 
 */

#include "optim.h"

/* Common Block Declarations */

struct
{
  double t0, tf, dti, dtf, ermx;
  int iu[5], nuc, nuv, ilin, nti, ntf, ny, nea, itmx, nex, nob, ntob, ntobi,
    nitu, ndtu;
} icsez_;

#define icsez_1 icsez_

struct
{
  int nitv, nrtv, ndtv;
} nird_;

#define nird_1 nird_


int optim_icse (int *ind, int *nu, double *u, double *co, double *g, int *itv,
		double *rtv, double *dtv, U_fp icsef, U_fp icsec2, U_fp icsei)
{
  /* System generated locals */
  int i__1, i__2;
  /* Local variables */
  int lech, lcof, indi, lobs, ltob, ldmy, lyob, ldtu, litu, mdtv, mitv, lsmy,
    ldif1, ldif2, ldif3;
  int lipv1, lipv2, mdtv1, mdtv2, mitv1, mitv2, i__, ludep, litob, loldp,
    lyold, lytob, ldtvt, lyerr, lyint, litvt, lytot, lb, ld, lf, lp, ly,
    lsmold, loldmu, lp0, ly0, lob, ldm, lfu, lui, nui, lfy, lgt, lc2y, ly0u;

  /* Parameter adjustments */
  --g;
  --u;
  --itv;
  --rtv;
  --dtv;

  /* Function Body */
  if (icsez_1.iu[1] > 0)
    {
      /*Computing MIN 
       */
      i__1 = *nu, i__2 = icsez_1.nuc + 1;
      lui = Min (i__1, i__2);
    }
  if (icsez_1.iu[0] > 0)
    {
      lui = 1;
    }
  nui =
    icsez_1.iu[0] * icsez_1.nuc + icsez_1.iu[1] * icsez_1.nuv * (icsez_1.nti +
								 icsez_1.ntf +
								 1);
  /* 
   *    decoupage de itv 
   *    nitu longueur de itu tableau de travail entier reserve 
   *    a l'utilisateur 
   *    nitvt longueur de itvt tableau de travail entier de 
   *    icse1 et icse2 
   * 
   */
  litu = 1;
  litvt = litu + icsez_1.nitu;
  /* 
   *    decoupage de dtv 
   *    ndtu longueur de dtu tableau de travail double precision 
   *    reserve a l'utilisateur 
   *    ndtvt longueur de dtvt tableau de travail double precision 
   *    de icse1 et icse2 
   * 
   */
  ldtu = 1;
  ly0 = ldtu + icsez_1.ndtu;
  ltob = ly0 + icsez_1.ny;
  lobs = ltob + icsez_1.ntob;
  lob = lobs + icsez_1.nob * icsez_1.ny;
  lech = lob + icsez_1.nex * icsez_1.ntob * icsez_1.nob;
  lcof = lech + *nu;
  /*      ********************** Modif 88 
   */
  lb = lcof + icsez_1.nob * icsez_1.ntob;
  lfy = lb + icsez_1.ny;
  lfu = lfy + icsez_1.ny * icsez_1.ny;
  ludep = lfu + icsez_1.ny * (icsez_1.nuc + icsez_1.nuv);
  lytot = ludep + *nu;
  lf = lytot + icsez_1.ny * (icsez_1.nti + icsez_1.ntf);
  ldtvt = lf + icsez_1.ny;
  /* 
   *    decoupage de itvt pour icse1 
   * 
   */
  lipv1 = litvt;
  mitv1 = lipv1 + icsez_1.ny - 1;
  /* 
   *    decoupage de itvt pour icse2 
   * 
   */
  litob = litvt;
  lipv2 = litob + icsez_1.ntob;
  mitv2 = lipv2 + icsez_1.ny - 1;
  /* 
   */
  mitv = Max (mitv1, mitv2);
  /* 
   *    decoupage de dtvt pour icse1 
   * 
   */
  ldm = ldtvt;
  lyold = ldm + icsez_1.ny * icsez_1.ny;
  lsmold = lyold + icsez_1.ny;
  lyint = lsmold + icsez_1.ny;
  lyerr = lyint + icsez_1.ny;
  ldif1 = lyerr + icsez_1.ny;
  ldif2 = ldif1 + icsez_1.ny;
  ldif3 = ldif2 + icsez_1.ny;
  mdtv1 = ldif3 + icsez_1.ny - 1;
  /* 
   *    decoupage de dtvt pour icse2 
   * 
   */
  lytob = ldtvt;
  lc2y = lytob + icsez_1.ny * icsez_1.ntob;
  ly0u = lc2y + icsez_1.ny * icsez_1.ntob;
  ldmy = ly0u + icsez_1.ny * *nu;
  lsmy = ldmy + icsez_1.ny * icsez_1.ny;
  loldmu = lsmy + icsez_1.ny * icsez_1.ny;
  ly = loldmu + icsez_1.ny * (icsez_1.nuc + icsez_1.nuv);
  loldp = ly + icsez_1.ny;
  lp = loldp + icsez_1.ny;
  lp0 = lp + icsez_1.ny;
  lgt = lp0 + icsez_1.ny;
  /*Computing MAX 
   */
  i__1 = icsez_1.nuc + icsez_1.nuv;
  lyob = lgt + Max (i__1, nui);
  ld = lyob + icsez_1.nob * icsez_1.ntob;
  mdtv2 = ld + icsez_1.nob - 1;
  /* 
   */
  mdtv = Max (mdtv1, mdtv2);
  if (mitv > nird_1.nitv || mdtv > nird_1.ndtv)
    {
      if (nird_1.nitv + nird_1.ndtv > 0)
	{
	  Sciprintf("icse: taille des tableaux itv,dtv insuffisante\n");
	  Sciprintf("\tvaleurs minimales %d et %d\n",mitv,mdtv);
	}
      nird_1.nitv = mitv;
      nird_1.ndtv = mdtv;
      return 0;
    }
  i__1 = *nu;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      dtv[ludep + i__ - 1] = u[i__];
      u[i__] = dtv[lech + i__ - 1] * u[i__];
      /* L10: */
    }
  /* 
   *    etat initial dependant du controle 
   * 
   */
  if (icsez_1.iu[0] > 0)
    {
      indi = 1;
      (*icsei) (&indi, &nui, &u[lui], &dtv[ly0], &dtv[ly0u], &itv[litu],
		&dtv[ldtu], &icsez_1.t0, &icsez_1.tf, &icsez_1.dti,
		&icsez_1.dtf, &icsez_1.ermx, icsez_1.iu, &icsez_1.nuc,
		&icsez_1.nuv, &icsez_1.ilin, &icsez_1.nti, &icsez_1.ntf,
		&icsez_1.ny, &icsez_1.nea, &icsez_1.itmx, &icsez_1.nex,
		&icsez_1.nob, &icsez_1.ntob, &icsez_1.ntobi, &icsez_1.nitu,
		&icsez_1.ndtu);
      if (indi <= 0)
	{
	  *ind = indi;
	  return 0;
	}
    }
  /* 
   *appel de icse1 
   *but 
   *    icse1 resout les systemes dynamiques du type: 
   *    dtu        double precision (ndtu) 
   *               tableau de travail double precision reserve 
   *               a l'utilisateur 
   * 
   *    enfin: 
   * 
   *    kt         entier 
   *               indice de comptage des pas de temps 
   * 
   *    dt         double precision 
   *               pas de temps,egal a dti ou a dtf 
   * 
   *    dtinv      double precision 
   *               dtinv=1/dt 
   * 
   *    t          double precision 
   *               instant(a l'instant t on travaille sur [t-dt,t]) 
   * 
   *    told       double precision 
   *               instant anterieur a t:told=t-dt 
   * 
   *    indf       entier 
   *               indicateur figurant dans la liste d'appel de icsef 
   * 
   *    it         entier 
   *               indice de comptage des corrections 
   * 
   *    err        double precision 
   *               norme l2 de dif2 
   * 
   */
  optim_icse1 (ind, nu, &u[1], (U_fp) icsef, &dtv[ly0], &dtv[lytot], &dtv[lf],
	       &dtv[lb], &dtv[lfy], &dtv[lfu], &itv[lipv1], &dtv[ldm],
	       &dtv[lyold], &dtv[lsmold], &dtv[lyint], &dtv[lyerr],
	       &dtv[ldif1], &dtv[ldif2], &dtv[ldif3], &itv[litu], &dtv[ldtu],
	       &icsez_1.t0, &icsez_1.tf, &icsez_1.dti, &icsez_1.dtf,
	       &icsez_1.ermx, icsez_1.iu, &icsez_1.nuc, &icsez_1.nuv,
	       &icsez_1.ilin, &icsez_1.nti, &icsez_1.ntf, &icsez_1.ny,
	       &icsez_1.nea, &icsez_1.itmx, &icsez_1.nex, &icsez_1.nob,
	       &icsez_1.ntob, &icsez_1.ntobi, &icsez_1.nitu, &icsez_1.ndtu);
  /* 
   */
  if (*ind <= 0)
    {
      return 0;
    }
  /* 
   *               stockage de df/du a l'instant posterieur 
   * 
   *    y          double precision (ny) 
   *               stockage de la valeur calculee de l'etat a un pas de 
   *               temps 
   * 
   *    oldp       double precision (ny) 
   *               stockage de la valeur calculee de l'etat adjoint au 
   *               pas de temps posterieur 
   * 
   *    p          double precision (ny) 
   *               stockage de la valeur calculee de l'etat adjoint a 
   *               un pas de temps 
   * 
   *    p0         double precision (ny) 
   *               etape dans le calcul des seconds membres des 
   *               systemes lineaires donnant l'etat adjoint dicretise 
   * 
   *    gt         double precision (nu) 
   *               stockage de la contribution au gradient a chaque 
   *               pas de temps 
   * 
   *    yob,d      figurent dans la liste d'appel de icsec2,qui calcule 
   *               le cout quadratique dans le cas d'un observateur 
   *               lineaire 
   * 
   *    itu        entier (nitu) 
   *               tableau de travail entier reserve a l'utlisateur 
   * 
   *    dtu        double precision (ndtu) 
   *               tableau de travail double precision reserve 
   *               a l'utilisateur 
   * 
   *    enfin: 
   * 
   *    kt         entier 
   *               indice de comptage des pas de temps 
   * 
   *    ktob       entier 
   *               indice de comptage des instants de mesure 
   * 
   *    dt         double precision 
   *               pas de temps,egal a dti ou a dtf 
   * 
   *    dt2        double precision 
   *               dt2=dt/2 
   * 
   *    dt2new      double precision 
   *               dt2 a l'instant posterieur 
   * 
   *    t          double precision 
   *               instant (a l'instant t on travaille sur [t,t+dt]) 
   * 
   *    c2         double precision 
   *               stockage d'un calcul inutile du cout 
   *               (au cas ou l'on n'utiliserait pas indc) 
   * 
   *    indf       entier 
   *               indicateur figurant dans la liste d'appel de icsef 
   * 
   *    indi       entier 
   *               indicateur figurant dans la liste d'appel de icsei 
   * 
   *    nui        entier 
   *               nombre de parametres dont depend l'etat initial 
   *               (figure dans la liste d'appel de icsei) 
   * 
   */
  optim_icse2 (ind, nu, &u[1], co, &g[1], (U_fp) icsef, (U_fp) icsec2,
	       (U_fp) icsei, &dtv[ly0], &dtv[ltob], &dtv[lobs], &dtv[lob],
	       &dtv[lytot], &dtv[lf], &dtv[lb], &dtv[lfy], &dtv[lfu],
	       &itv[lipv2], &itv[litob], &dtv[lcof], &dtv[lytob], &dtv[lc2y],
	       &dtv[ly0u], &dtv[ldmy], &dtv[lsmy], &dtv[loldmu], &dtv[ly],
	       &dtv[loldp], &dtv[lp], &dtv[lp0], &dtv[lgt], &dtv[lyob],
	       &dtv[ld], &itv[litu], &dtv[ldtu], &icsez_1.t0, &icsez_1.tf,
	       &icsez_1.dti, &icsez_1.dtf, &icsez_1.ermx, icsez_1.iu,
	       &icsez_1.nuc, &icsez_1.nuv, &icsez_1.ilin, &icsez_1.nti,
	       &icsez_1.ntf, &icsez_1.ny, &icsez_1.nea, &icsez_1.itmx,
	       &icsez_1.nex, &icsez_1.nob, &icsez_1.ntob, &icsez_1.ntobi,
	       &icsez_1.nitu, &icsez_1.ndtu,"","","");
  i__1 = *nu;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      g[i__] = dtv[lech + i__ - 1] * g[i__];
      u[i__] = dtv[ludep + i__ - 1];
      /* L20: */
    }
  return 0;
}



int optim_icse2 (int *ind, int *nu, double *u, double *co, double *g, U_fp icsef,
		 U_fp icsec2, U_fp icsei, double *y0, double *tob, double *obs,
		 double *ob, double *ytot, double *f, double *b, double *fy,
		 double *fu, int *ipv2, int *itob, double *cof, double *ytob,
		 double *c2y, double *y0u, double *dmy, double *smy,
		 double *oldmu, double *y, double *oldp, double *p, double *p0,
		 double *gt, double *yob, double *d__, int *itu, double *dtu,
		 double *t0, double *tf, double *dti, double *dtf, double *ermx,
		 int *iu, int *nuc, int *nuv, int *ilin, int *nti, int *ntf,
		 int *ny, int *nea, int *itmx, int *nex, int *nob, int *ntob,
		 int *ntobi, int *nitu, int *ndtu, char *nomf, char *nomc,
		 char *nomi)
{
  double c_b2 = 0.;
  int c__1 = 1;

  /* System generated locals */
  int obs_dim1, obs_offset, ob_dim1, ob_dim2, ob_offset, ytot_dim1,
    ytot_offset, fy_dim1, fy_offset, fu_dim1, fu_offset, cof_dim1, cof_offset,
    ytob_dim1, ytob_offset, c2y_dim1, c2y_offset, y0u_dim1, y0u_offset,
    dmy_dim1, dmy_offset, smy_dim1, smy_offset, oldmu_dim1, oldmu_offset,
    yob_dim1, yob_offset, i__1, i__2;

  /* Local variables */
  int indc, indf, indi;
  int ktob, info;
  int i__, j;
  double t;
  double dt2new, dt;
  int kt;
  double dt2;
  int lui, nui, luv;

  /* 
   *    sous programme de icse.fortran:calcul du gradient par 
   *    integration du systeme adjoint 
   * 
   * 
   * 
   *    Initialisation 
   * 
   */
  /* Parameter adjustments */
  --gt;
  --g;
  --u;
  --iu;
  --p0;
  --p;
  --oldp;
  --y;
  oldmu_dim1 = *ny;
  oldmu_offset = oldmu_dim1 + 1;
  oldmu -= oldmu_offset;
  smy_dim1 = *ny;
  smy_offset = smy_dim1 + 1;
  smy -= smy_offset;
  dmy_dim1 = *ny;
  dmy_offset = dmy_dim1 + 1;
  dmy -= dmy_offset;
  y0u_dim1 = *ny;
  y0u_offset = y0u_dim1 + 1;
  y0u -= y0u_offset;
  --ipv2;
  fu_dim1 = *ny;
  fu_offset = fu_dim1 + 1;
  fu -= fu_offset;
  fy_dim1 = *ny;
  fy_offset = fy_dim1 + 1;
  fy -= fy_offset;
  --b;
  --f;
  ytot_dim1 = *ny;
  ytot_offset = ytot_dim1 + 1;
  ytot -= ytot_offset;
  --y0;
  --d__;
  obs_dim1 = *nob;
  obs_offset = obs_dim1 + 1;
  obs -= obs_offset;
  yob_dim1 = *nob;
  yob_offset = yob_dim1 + 1;
  yob -= yob_offset;
  c2y_dim1 = *ny;
  c2y_offset = c2y_dim1 + 1;
  c2y -= c2y_offset;
  ytob_dim1 = *ny;
  ytob_offset = ytob_dim1 + 1;
  ytob -= ytob_offset;
  cof_dim1 = *nob;
  cof_offset = cof_dim1 + 1;
  cof -= cof_offset;
  --itob;
  ob_dim1 = *nex;
  ob_dim2 = *ntob;
  ob_offset = ob_dim1 * (ob_dim2 + 1) + 1;
  ob -= ob_offset;
  --tob;
  --itu;
  --dtu;

  /* Function Body */
  nsp_dset (nu, &c_b2, &g[1], &c__1);
  nsp_dset (ny, &c_b2, &p[1], &c__1);
  kt = *nti + *ntf;
  ktob = *ntob;
  /*    lui et nui servent quand l'etat initial depend du controle 
   */
  if (iu[2] > 0)
    {
      /*Computing MIN 
       */
      i__1 = *nu, i__2 = *nuc + 1;
      lui = Min (i__1, i__2);
    }
  if (iu[1] > 0)
    {
      lui = 1;
    }
  nui = iu[1] * *nuc + iu[2] * *nuv * (*nti + *ntf + 1);
  /* 
   * 
   *    Calcul de itob,vecteur des indices des instants de mesure 
   *    a partir de tob,vecteur des instants de mesure 
   *    itob(j) est l'entier le plus proche de tob(j)/dt 
   * 
   */
  i__1 = *ntobi;
  for (j = 1; j <= i__1; ++j)
    {
      /* L1: */
      itob[j] = (int) ((tob[j] - *t0) / *dti + .5);
    }
  if (*ntobi < *ntob)
    {
      itob[*ntobi + 1] =
	*nti + (int) ((tob[*ntobi + 1] - *t0 - *nti * *dti) / *dtf + .5);
    }
  if (*ntobi + 1 < *ntob)
    {
      i__1 = *ntob;
      for (j = *ntobi + 2; j <= i__1; ++j)
	{
	  /* L2: */
	  itob[j] =
	    itob[*ntobi + 1] + (int) ((tob[j] - tob[*ntobi + 1]) / *dtf + .5);
	}
    }
  /* 
   *    Ecriture de ytob tableau des valeurs de l'etat 
   *    aux instants de mesure 
   * 
   */
  i__1 = *ntob;
  for (j = 1; j <= i__1; ++j)
    {
      i__2 = *ny;
      for (i__ = 1; i__ <= i__2; ++i__)
	{
	  /* L10: */
	  C2F(dcopy) (ny, &ytot[itob[j] * ytot_dim1 + 1], &c__1,
		      &ytob[j * ytob_dim1 + 1], &c__1);
	}
    }
  /* 
   *    Si ind=2,on calcule seulement le cout 
   *    Si ind=3,on calcule seulement le gradient 
   *    Si ind=4,on calcule le cout et le gradient 
   * 
   */
  if (*ind != 3)
    {
      indc = 1;
      (*icsec2) (&indc, nu, &tob[1], &obs[obs_offset], &cof[cof_offset],
		 &ytob[ytob_offset], &ob[ob_offset], &u[1], co,
		 &c2y[c2y_offset], &g[1], &yob[yob_offset], &d__[1], &itu[1],
		 &dtu[1], t0, tf, dti, dtf, ermx, &iu[1], nuc, nuv, ilin, nti,
		 ntf, ny, nea, itmx, nex, nob, ntob, ntobi, nitu, ndtu, nomf,
		 nomc, nomi, 6L, 6L, 6L);
      if (indc <= 0)
	{
	  *ind = indc;
	  return 0;
	}
    }
  if (*ind == 2)
    {
      return 0;
    }
  /* 
   *    Calcul du gradient du cout en utilisant l'etat adjoint 
   *    discretise: 
   *    calcul de la derivee partielle c2y du cout par rapport 
   *    a l'etat 
   * 
   */
  indc = 2;
  (*icsec2) (&indc, nu, &tob[1], &obs[obs_offset], &cof[cof_offset],
	     &ytob[ytob_offset], &ob[ob_offset], &u[1], co, &c2y[c2y_offset],
	     &g[1], &yob[yob_offset], &d__[1], &itu[1], &dtu[1], t0, tf, dti,
	     dtf, ermx, &iu[1], nuc, nuv, ilin, nti, ntf, ny, nea, itmx, nex,
	     nob, ntob, ntobi, nitu, ndtu, nomf, nomc, nomi, 6L, 6L, 6L);
  if (indc <= 0)
    {
      *ind = indc;
      return 0;
    }
  /* 
   *    +Evaluations successives de la contribution au gradient 
   *    a chaque pas de temps a l'aide de l'etat adjoint(non stocke) 
   * 
   */
  for (kt = *nti + *ntf; kt >= 1; --kt)
    {
      /* 
       *    *Calcul de l'etat adjoint p_kt au pas kt 
       *      Calcul du second membre p 
       *      Initialisation: 
       *      nt=nti+ntf 
       *      si kt=nt,p est nul 
       *      si kt<nt,on a au depart p=p_kt+1,etat adjoint au pas kt+1; 
       *      on prend p=p0,ou p0=(smy)t*I*p avec smy=Id+(dt/2).dfy(t,y,u) 
       *      avec dt=optim_t(kt+1)-t_kt,y=y_kt,t=t_kt,I designant la matrice 
       *      diagonale d'ordre ny dont les nea premiers coefficients 
       *      valent 0 et les autres 1 et Id designant la matrice 
       *      identite d'ordre ny; 
       *      si le systeme est affine avec partie lineaire autonome 
       *      (ilin=2) smy est calculee seulement aux premiers pas de temps 
       *      pour dti et dtf,sinon (ilin=0 ou 1) smy est calculee a chaque 
       *      pas de temps 
       * 
       *      stockage de l'etat adjoint et de dt/2 au pas kt+1 
       * 
       */
      C2F(dcopy) (ny, &p[1], &c__1, &oldp[1], &c__1);
      /*Computing MIN 
       */
      i__2 = *nu, i__1 = *nuc + 1 + kt * *nuv;
      luv = Min (i__2, i__1);
      /* 
       *      calcul de y=y_kt,dt=optim_t(kt+1)-t_kt,dt2=dt/2, 
       *                dt2new=(t_kt-optim_t(kt-1))/2 
       * 
       */
      C2F(dcopy) (ny, &ytot[kt * ytot_dim1 + 1], &c__1, &y[1], &c__1);
      /* 
       */
      if (kt < *nti)
	{
	  t = kt * *dti + *t0;
	  dt = *dti;
	}
      else
	{
	  t = *nti * *dti + (kt - *nti) * *dtf + *t0;
	  dt = *dtf;
	}
      dt2 = dt / 2.;
      if (kt != *nti)
	{
	  dt2new = dt2;
	}
      else
	{
	  dt2new = *dti / 2.;
	}
      /* 
       *      Dans le cas ilin<=1, 
       *        calcul de fy=dfy(t,y,u) puis de smy=Id+(dt/2).dmy 
       *        lorsque kt<(nti+ntf) 
       *      Sinon (ilin>1),fy=dfy a ete calcule dans icse1 
       * 
       * 
       */
      if (*ilin <= 1)
	{
	  indf = 2;
	  (*icsef) (&indf, &t, &y[1], &u[1], &u[luv], &f[1], &fy[fy_offset],
		    &fu[fu_offset], &b[1], &itu[1], &dtu[1], t0, tf, dti, dtf,
		    ermx, &iu[1], nuc, nuv, ilin, nti, ntf, ny, nea, itmx,
		    nex, nob, ntob, ntobi, nitu, ndtu, nomf, nomc, nomi, 6L,
		    6L, 6L);
	  if (indf <= 0)
	    {
	      *ind = indf;
	      return 0;
	    }
	}
      /* 
       */
      if (kt != *nti + *ntf)
	{
	  if (*ilin <= 1 || kt == *nti + *ntf - 1 || kt == *nti - 1)
	    {
	      i__2 = *ny;
	      for (i__ = 1; i__ <= i__2; ++i__)
		{
		  i__1 = *ny;
		  for (j = 1; j <= i__1; ++j)
		    {
		      /* L30: */
		      smy[i__ + j * smy_dim1] = dt2 * fy[i__ + j * fy_dim1];
		    }
		}
	      i__1 = *ny;
	      for (i__ = 1; i__ <= i__1; ++i__)
		{
		  /* L35: */
		  smy[i__ + i__ * smy_dim1] += 1.;
		}
	    }
	  /* 
	   *        calcul de p0=(smy)t*I*p puis p=p0 
	   * 
	   */
	  if (*nea > 0)
	    {
	      i__1 = *nea;
	      for (i__ = 1; i__ <= i__1; ++i__)
		{
		  /* L40: */
		  p[i__] = 0.;
		}
	    }
	  dmmul_scicos (&p[1], &c__1, &smy[smy_offset], ny, &p0[1], &c__1,
			&c__1, ny, ny);
	  /* 
	   */
	  C2F(dcopy) (ny, &p0[1], &c__1, &p[1], &c__1);
	}
      /* 
       *      Fin du calcul du second membre p 
       *        si kt=itob(ktob),on ajoute c2y(.,ktob) au second membre p 
       * 
       */
      if (ktob > 0)
	{
	  if (kt == itob[ktob])
	    {
	      i__1 = *ny;
	      for (i__ = 1; i__ <= i__1; ++i__)
		{
		  /* L50: */
		  p[i__] += c2y[i__ + ktob * c2y_dim1];
		}
	      --ktob;
	    }
	}
      /* 
       *      Calcul et factorisation de la matrice dmy du systeme 
       *      de l'etat adjoint 
       *      dmy=I-dt2new.dfy(t,y,u) avec dt2new=(t_kt-optim_t(kt-1))/2, 
       *      y=y_kt,t=t_kt,Idesignant la matrice diagonale d'ordre 
       *      ny dont les nea premiers coefficients valent 0 et les 
       *      autres 1; 
       *      si le systeme est affine avec partie lineaire autonome 
       *      (ilin=2) dmy est calculee et factorisee aux premiers 
       *      pas de temps pour dti et dtf,sinon (ilin=0 ou 1) dmy est 
       *      calculee et factorisee a chaque pas de temps 
       * 
       */
      if (*ilin <= 1 || kt == *nti + *ntf || kt == *nti)
	{
	  i__1 = *ny;
	  for (i__ = 1; i__ <= i__1; ++i__)
	    {
	      i__2 = *ny;
	      for (j = 1; j <= i__2; ++j)
		{
		  /* L60: */
		  dmy[i__ + j * dmy_dim1] = -dt2new * fy[i__ + j * fy_dim1];
		}
	    }
	  i__2 = *ny;
	  for (i__ = *nea + 1; i__ <= i__2; ++i__)
	    {
	      /* L65: */
	      dmy[i__ + i__ * dmy_dim1] += 1.;
	    }
	  C2F(dgefa) (&dmy[dmy_offset], ny, ny, &ipv2[1], &info);
	}
      /* 
       *      Resolution de (dmy)t*X=p,la solution s'appelant p 
       *      p est alors p_kt,etat adjoint au pas kt 
       * 
       */
      C2F(dgesl) (&dmy[dmy_offset], ny, ny, &ipv2[1], &p[1], &c__1);
      /* 
       *    *Calcul du gradient g au pas kt+1 
       *      calcul de la contribution gt au gradient au pas kt+1: 
       *      gt=(dt/2).(I*dfu(t_kt,y_kt,u)+dfu(t_kt+1,y_kt+1,u))t*p_kt+1 
       *      avec dt=optim_t(kt+1)-t_kt,I designant la matrice diagonale 
       *      d'ordre ny dont les nea premiers coefficients valent 0 
       *      et les autres 1 
       * 
       */
      if (*nuv > 0 || iu[3] == 1)
	{
	  indf = 3;
	  (*icsef) (&indf, &t, &y[1], &u[1], &u[luv], &f[1], &fy[fy_offset],
		    &fu[fu_offset], &b[1], &itu[1], &dtu[1], t0, tf, dti, dtf,
		    ermx, &iu[1], nuc, nuv, ilin, nti, ntf, ny, nea, itmx,
		    nex, nob, ntob, ntobi, nitu, ndtu, nomf, nomc, nomi, 6L,
		    6L, 6L);
	  if (indf <= 0)
	    {
	      *ind = indf;
	      return 0;
	    }
	  if (kt < *nti + *ntf)
	    {
	      i__2 = *nuc + *nuv;
	      dmmul_scicos (&oldp[1], &c__1, &oldmu[oldmu_offset], ny, &gt[1],
			    &c__1, &c__1, ny, &i__2);
	      i__2 = *nuc + *nuv;
	      C2F(dscal) (&i__2, &dt2, &gt[1], &c__1);
	      /*          le gradient g est la somme des contributions 
	       */
	      if (iu[3] > 0)
		{
		  nsp_dadd (*nuc, &gt[1], c__1, &g[1], c__1);
		}
	      if (*nuv > 0)
		{
		  /*Computing MIN 
		   */
		  i__2 = *nu, i__1 = *nuc + 1 + (kt + 1) * *nuv;
		  luv = Min (i__2, i__1);
		  nsp_dadd (*nuv, &gt[*nuc + 1], c__1, &g[luv], c__1);
		}
	      if (*nea > 0)
		{
		  i__2 = *nea;
		  for (i__ = 1; i__ <= i__2; ++i__)
		    {
		      /* L70: */
		      oldp[i__] = 0.;
		    }
		}
	      i__2 = *nuc + *nuv;
	      dmmul_scicos (&oldp[1], &c__1, &fu[fu_offset], ny, &gt[1], &c__1,
			    &c__1, ny, &i__2);
	      i__2 = *nuc + *nuv;
	      C2F(dscal) (&i__2, &dt2, &gt[1], &c__1);
	      /*          le gradient g est la somme des contributions 
	       */
	      if (iu[3] > 0)
		{
		  nsp_dadd (*nuc, &gt[1], c__1, &g[1], c__1);
		}
	      if (*nuv > 0)
		{
		  /*Computing MIN 
		   */
		  i__2 = *nu, i__1 = *nuc + 1 + kt * *nuv;
		  luv = Min (i__2, i__1);
		  nsp_dadd (*nuv, &gt[*nuc + 1], c__1, &g[luv], c__1);
		}
	    }
	  /* 
	   *        stockage de dfu(t_kt,y_kt,u) dans oldmu 
	   * 
	   */
	  i__2 = *ny * (*nuc + *nuv);
	  C2F(dcopy) (&i__2, &fu[fu_offset], &c__1, &oldmu[oldmu_offset],
		      &c__1);
	  /* 
	   *    *On passe au pas de temps suivant:kt-1,sauf si kt=1,auquel cas 
	   *      on calcule la contribution gt au gradient au pas kt=1 et 
	   *      on l'ajoute a g;on a: 
	   *      gt=(dt/2).(I*dfu(t0,y0,u)+dfu(t_1,y_1,u))t*p_1,avec 
	   *      dt=t_1-t0=dti car nti n'est jamais nul par convention 
	   * 
	   */
	  if (kt == 1)
	    {
	      t = *t0;
	      dt2 = *dti / 2.;
	      C2F(dcopy) (ny, &y0[1], &c__1, &y[1], &c__1);
	      indf = 3;
	      (*icsef) (&indf, &t, &y[1], &u[1], &u[luv], &f[1],
			&fy[fy_offset], &fu[fu_offset], &b[1], &itu[1],
			&dtu[1], t0, tf, dti, dtf, ermx, &iu[1], nuc, nuv,
			ilin, nti, ntf, ny, nea, itmx, nex, nob, ntob, ntobi,
			nitu, ndtu, nomf, nomc, nomi, 6L, 6L, 6L);
	      if (indf <= 0)
		{
		  *ind = indf;
		  return 0;
		}
	      i__2 = *nuc + *nuv;
	      dmmul_scicos (&p[1], &c__1, &oldmu[oldmu_offset], ny, &gt[1],
			    &c__1, &c__1, ny, &i__2);
	      i__2 = *nuc + *nuv;
	      C2F(dscal) (&i__2, &dt2, &gt[1], &c__1);
	      /*          le gradient g est la somme des contributions 
	       */
	      if (iu[3] > 0)
		{
		  nsp_dadd (*nuc, &gt[1], c__1, &g[1], c__1);
		}
	      if (*nuv > 0)
		{
		  /*Computing MIN 
		   */
		  i__2 = *nu, i__1 = *nuc + 1 + *nuv;
		  luv = Min (i__2, i__1);
		  nsp_dadd (*nuv, &gt[*nuc + 1], c__1, &g[luv], c__1);
		}
	      if (*nea > 0)
		{
		  i__2 = *nea;
		  for (i__ = 1; i__ <= i__2; ++i__)
		    {
		      /* L90: */
		      p[i__] = 0.;
		    }
		}
	      i__2 = *nuc + *nuv;
	      dmmul_scicos (&p[1], &c__1, &fu[fu_offset], ny, &gt[1], &c__1,
			    &c__1, ny, &i__2);
	      i__2 = *nuc + *nuv;
	      C2F(dscal) (&i__2, &dt2, &gt[1], &c__1);
	      /*          le gradient g est la somme des contributions 
	       */
	      if (iu[3] > 0)
		{
		  nsp_dadd (*nuc, &gt[1], c__1, &g[1], c__1);
		}
	      if (*nuv > 0)
		{
		  /*Computing MIN 
		   */
		  i__2 = *nu, i__1 = *nuc + 1;
		  luv = Min (i__2, i__1);
		  nsp_dadd (*nuv, &gt[*nuc + 1], c__1, &g[luv], c__1);
		}
	    }
	}
      /* L100: */
    }
  /* 
   *    gradient par rapport au controle initial 
   * 
   */
  if (Max (iu[1], iu[2]) > 0)
    {
      /* 
       *    calcul de l'etat adjoint initial p0 
       * 
       */
      indf = 2;
      (*icsef) (&indf, &t, &y[1], &u[1], &u[luv], &f[1], &fy[fy_offset],
		&fu[fu_offset], &b[1], &itu[1], &dtu[1], t0, tf, dti, dtf,
		ermx, &iu[1], nuc, nuv, ilin, nti, ntf, ny, nea, itmx, nex,
		nob, ntob, ntobi, nitu, ndtu, nomf, nomc, nomi, 6L, 6L, 6L);
      if (indf == 0)
	{
	  *ind = indf;
	  return 0;
	}
      i__2 = *ny;
      for (i__ = 1; i__ <= i__2; ++i__)
	{
	  i__1 = *ny;
	  for (j = 1; j <= i__1; ++j)
	    {
	      /* L120: */
	      smy[i__ + j * smy_dim1] = dt2 * fy[i__ + j * fy_dim1];
	    }
	}
      i__1 = *ny;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  /* L125: */
	  smy[i__ + i__ * smy_dim1] += 1.;
	}
      if (*nea > 0)
	{
	  i__1 = *nea;
	  for (i__ = 1; i__ <= i__1; ++i__)
	    {
	      /* L130: */
	      p[i__] = 0.;
	    }
	}
      dmmul_scicos (&p[1], &c__1, &smy[smy_offset], ny, &p0[1], &c__1, &c__1,
		    ny, ny);
      /*    incrementation du gradient 
       */
      indi = 2;
      (*icsei) (&indi, &nui, &u[lui], &y0[1], &y0u[y0u_offset], &itu[1],
		&dtu[1], t0, tf, dti, dtf, ermx, &iu[1], nuc, nuv, ilin, nti,
		ntf, ny, nea, itmx, nex, nob, ntob, ntobi, nitu, ndtu, nomf,
		nomc, nomi, 6L, 6L, 6L);
      if (indi <= 0)
	{
	  *ind = indi;
	  return 0;
	}
      dmmul_scicos (&p0[1], &c__1, &y0u[y0u_offset], ny, &gt[1], &c__1, &c__1,
		    &nui, &nui);
      i__1 = nui;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  /* L150: */
	  g[lui + i__ - 1] += gt[i__];
	}
      /* 
       */
    }
  return 0;
}


/* 
 *    sous programme de icse.fortran:calcul de l'etat 
 *    Copyright INRIA 
 */

int optim_icse1 (int *ind, int *nu, double *u, U_fp icsef, double *y0,
		 double *ytot, double *f, double *b, double *fy, double *fu,
		 int *ipv1, double *dm, double *yold, double *smold, double *yint,
		 double *yerr, double *dif1, double *dif2, double *dif3, int *itu,
		 double *dtu, double *t0, double *tf, double *dti, double *dtf,
		 double *ermx, int *iu, int *nuc, int *nuv, int *ilin, int *nti,
		 int *ntf, int *ny, int *nea, int *itmx, int *nex, int *nob,
		 int *ntob, int *ntobi, int *nitu, int *ndtu)
{
  int c__1 = 1;
  double c_b17 = .5;
  int c__0 = 0;
  /* System generated locals */
  int ytot_dim1, ytot_offset, fy_dim1, fy_offset, fu_dim1, fu_offset, dm_dim1,
    dm_offset, i__1, i__2, i__3;
  double d__1;

  /* Local variables */
  int indf, info;
  double told;
  int i__, j;
  double t;
  double dtinv;
  double dt;
  int it, kt;
  double err;
  int luv;


  /* Parameter adjustments */
  --u;
  --iu;
  --dif3;
  --dif2;
  --dif1;
  --yerr;
  --yint;
  --smold;
  --yold;
  dm_dim1 = *ny;
  dm_offset = dm_dim1 + 1;
  dm -= dm_offset;
  --ipv1;
  fu_dim1 = *ny;
  fu_offset = fu_dim1 + 1;
  fu -= fu_offset;
  fy_dim1 = *ny;
  fy_offset = fy_dim1 + 1;
  fy -= fy_offset;
  --b;
  --f;
  ytot_dim1 = *ny;
  ytot_offset = ytot_dim1 + 1;
  ytot -= ytot_offset;
  --y0;
  --itu;
  --dtu;

  /* Function Body */
  t = *t0;
  C2F(dcopy) (ny, &y0[1], &c__1, &yold[1], &c__1);
  /* 
   *    Resolutions successives du systeme d'etat discretise a chaque 
   *    pas de temps 
   * 
   */
  i__1 = *nti + *ntf;
  for (kt = 1; kt <= i__1; ++kt)
    {
      /* 
       *    *Calcul,puis factorisation de la matrice dm du systeme: 
       *    on a au depart yold=y_kt-1,etat au pas kt-1;alors 
       *    dm=(1/dt).I-(1/2).dfy(told,yold,u) avec dt=t_kt-optim_t(kt-1), 
       *    told=optim_t(kt-1),I designant la matrice diagonale d'ordre ny 
       *    dont les nea premiers coefficients diagonaux valent 0 et les 
       *    autres 1 
       *    si le systeme est affine avec partie lineaire autonome (ilin=2) 
       *    dm est calculee et factorisee seulement aux premiers pas de 
       *    temps pour dti et dtf,sinon (ilin=0 ou 1) dm est calculee et 
       *    factorisee a chaque pas de temps 
       * 
       *    stockage de l'instant au pas kt-1 
       * 
       *Computing MIN 
       */
      i__2 = *nu, i__3 = *nuc + 1 + (kt - 1) * *nuv;
      luv = Min (i__2, i__3);
      told = t;
      /* 
       *    calcul de t=t_kt,dt=t_kt-optim_t(kt-1),dtinv=1/dt 
       * 
       */
      if (kt <= *nti)
	{
	  t = kt * *dti + *t0;
	  dt = *dti;
	}
      else
	{
	  t = *nti * *dti + (kt - *nti) * *dtf + *t0;
	  dt = *dtf;
	}
      dtinv = 1. / dt;
      /* 
       *    calcul et factorisation de dm=dtinv.I-(1/2).dfy(told,yold,u) 
       *    I designant la matrice diagonale d'ordre ny dont les nea 
       *    premiers coefficients diagonaux valent 0 et les autres 1 
       * 
       *    fy=dfy(told,yold,u) n'est calcule que pour kt=1 lorsque 
       *    ilin>1 
       *    dm=dtinv.I-(1/2).fy n'est calcule que pour kt=1 ou kt=nti+1 
       *    lorsque ilin>1 
       * 
       */
      if (kt == 1 || kt == *nti + 1 || *ilin <= 1)
	{
	  indf = 2;
	  if (kt == 1 || *ilin <= 1)
	    {
	      (*icsef) (&indf, &told, &yold[1], &u[1], &u[luv], &f[1],
			&fy[fy_offset], &fu[fu_offset], &b[1], &itu[1],
			&dtu[1], t0, tf, dti, dtf, ermx, &iu[1], nuc, nuv,
			ilin, nti, ntf, ny, nea, itmx, nex, nob, ntob, ntobi,
			nitu, ndtu);
	    }
	  if (indf <= 0)
	    {
	      *ind = indf;
	      return 0;
	    }
	  i__2 = *ny;
	  for (i__ = 1; i__ <= i__2; ++i__)
	    {
	      i__3 = *ny;
	      for (j = 1; j <= i__3; ++j)
		{
		  /* L10: */
		  dm[i__ + j * dm_dim1] = -fy[i__ + j * fy_dim1] / 2.;
		}
	    }
	  i__3 = *ny;
	  for (i__ = *nea + 1; i__ <= i__3; ++i__)
	    {
	      /* L20: */
	      dm[i__ + i__ * dm_dim1] += dtinv;
	    }
	  C2F(dgefa) (&dm[dm_offset], ny, ny, &ipv1[1], &info);
	}
      /* 
       *    *Calcul de l'etat y_kt au pas kt: 
       *      Initialisation du nombre d'iterations dans l'algorithme 
       * 
       */
      it = 1;
      /* 
       *      Prediction du deplacement:dif1 (Euler explicite) 
       *      dif1=dt.(I*f(told,yold,u)) 
       *      et stockage de I*f(told,yold,u) dans smold,I designant la 
       *      matrice diagonale d'ordre ny dont les nea premiers 
       *      coefficients diagonaux valent 0 et les autres 1 
       * 
       *      si kt=1,smold=f(told,yold,u) 
       *      sinon,smold=f(told,yold,u) a ete calcule au pas kt-1 
       *      sous le nom dif1=f(t,yerr,u) 
       * 
       */
      if (kt == 1)
	{
	  indf = 1;
	  (*icsef) (&indf, &told, &yold[1], &u[1], &u[luv], &smold[1],
		    &fy[fy_offset], &fu[fu_offset], &b[1], &itu[1], &dtu[1],
		    t0, tf, dti, dtf, ermx, &iu[1], nuc, nuv, ilin, nti, ntf,
		    ny, nea, itmx, nex, nob, ntob, ntobi, nitu, ndtu);
	  if (indf <= 0)
	    {
	      *ind = indf;
	      return 0;
	    }
	}
      /* 
       *      smold=I*smold,puis dif1=dt.smold 
       * 
       */
      if (*nea > 0)
	{
	  i__3 = *nea;
	  for (i__ = 1; i__ <= i__3; ++i__)
	    {
	      /* L30: */
	      smold[i__] = 0.;
	    }
	}
      C2F(dcopy) (ny, &smold[1], &c__1, &dif1[1], &c__1);
      C2F(dscal) (ny, &dt, &dif1[1], &c__1);
      /* 
       *      Calcul de l'erreur:dif2=(1/2).(smold+f(t,yold+dif1,u))-dif1/dt 
       *        dif2=f(t,yint,u) ou yint=yold+dif1 
       * 
       *Computing MIN 
       */
      i__3 = *nu, i__2 = *nuc + 1 + kt * *nuv;
      luv = Min (i__3, i__2);
      C2F(dcopy) (ny, &yold[1], &c__1, &yint[1], &c__1);
      nsp_dadd (*ny, &dif1[1], c__1, &yint[1], c__1);
      indf = 1;
      (*icsef) (&indf, &t, &yint[1], &u[1], &u[luv], &dif2[1], &fy[fy_offset],
		&fu[fu_offset], &b[1], &itu[1], &dtu[1], t0, tf, dti, dtf,
		ermx, &iu[1], nuc, nuv, ilin, nti, ntf, ny, nea, itmx, nex,
		nob, ntob, ntobi, nitu, ndtu);
      if (indf <= 0)
	{
	  *ind = indf;
	  return 0;
	}
      /* 
       *        dif2=(1/2).(smold+dif2) 
       * 
       */
      nsp_dadd (*ny, &smold[1], c__1, &dif2[1], c__1);
      C2F(dscal) (ny, &c_b17, &dif2[1], &c__1);
      /* 
       *        dif2=dif2-dif1/dt 
       * 
       */
      d__1 = -dtinv;
      C2F(daxpy) (ny, &d__1, &dif1[1], &c__1, &dif2[1], &c__1);
      /* 
       *      Calcul du nouvel ecart:dif3 
       *        initialisation:dif3=dif1 
       * 
       */
      C2F(dcopy) (ny, &dif1[1], &c__1, &dif3[1], &c__1);
      /* 
       *        resolution de dm*X=dif2,la solution X s'appelant dif2 
       * 
       */
    L50:
      C2F(dgesl) (&dm[dm_offset], ny, ny, &ipv1[1], &dif2[1], &c__0);
      /* 
       *        dif3=dif3+dif2 est le nouvel ecart 
       * 
       */
      nsp_dadd (*ny, &dif2[1], c__1, &dif3[1], c__1);
      /* 
       *      Calcul de la norme err de l'erreur dif2 apres correction 
       *      dif2=(1/2).(smold+f(t,yold+dif3,u))-I*(dif3/dt) 
       *      I designant la matrice diagonale d'ordre ny dont les nea 
       *      premiers coefficients diagonaux valent 0 et les autres 1 
       * 
       *        dif1=f(t,yerr,u) ou yerr=yold+dif3 
       * 
       */
      C2F(dcopy) (ny, &yold[1], &c__1, &yerr[1], &c__1);
      nsp_dadd (*ny, &dif3[1], c__1, &yerr[1], c__1);
      /* 
       *        ermx<0 correspond au fait que l'utilisateur a choisi 
       *        de ne faire qu'une correction sans aucun test sur 
       *        l'erreur 
       */
      if (*ermx < 0.)
	{
	  goto L55;
	}
      /* 
       */
      indf = 1;
      (*icsef) (&indf, &t, &yerr[1], &u[1], &u[luv], &dif1[1], &fy[fy_offset],
		&fu[fu_offset], &b[1], &itu[1], &dtu[1], t0, tf, dti, dtf,
		ermx, &iu[1], nuc, nuv, ilin, nti, ntf, ny, nea, itmx, nex,
		nob, ntob, ntobi, nitu, ndtu);
      if (indf <= 0)
	{
	  *ind = indf;
	  return 0;
	}
      /* 
       *        dif2=dif1 
       * 
       */
      C2F(dcopy) (ny, &dif1[1], &c__1, &dif2[1], &c__1);
      /* 
       *        dif2=(1/2).(smold+dif2) 
       * 
       */
      nsp_dadd (*ny, &smold[1], c__1, &dif2[1], c__1);
      C2F(dscal) (ny, &c_b17, &dif2[1], &c__1);
      /* 
       *        dif2=dif2-I*(dif3/dt) 
       * 
       */
      i__3 = *ny - *nea;
      d__1 = -dtinv;
      C2F(daxpy) (&i__3, &d__1, &dif3[*nea + 1], &c__1, &dif2[*nea + 1],
		  &c__1);
      /* 
       *        err=norme l2 de dif2 
       * 
       */
      err = C2F(dnrm2) (ny, &dif2[1], &c__1);
      /* 
       *      Si err>ermx,on corrige a nouveau si possible,c'est a dire 
       *      si it<=itmx 
       * 
       */
      if (err > *ermx && *ilin == 0)
	{
	  ++it;
	  if (it > *itmx)
	    {
	      *ind = -1;
	      Sciprintf("icse: integration de l'etat impossible");
	      return 0;
	    }
	  goto L50;
	}
      /* 
       *      Si err<=ermx,yold=yerr est y_kt,etat au pas kt 
       *      (on avait calcule yerr=yold+dif3) 
       *      et ytot(.,kt)=yold 
       * 
       */
    L55:
      C2F(dcopy) (ny, &yerr[1], &c__1, &yold[1], &c__1);
      C2F(dcopy) (ny, &yold[1], &c__1, &ytot[kt * ytot_dim1 + 1], &c__1);
      /* 
       *      smold=dif1 
       *      (on avait calcule dif1=f(t,yerr,u) qui deviendra 
       *      f(told,yold,u) au pas kt+1) 
       * 
       */
      C2F(dcopy) (ny, &dif1[1], &c__1, &smold[1], &c__1);
      /* 
       *    *On passe au pas de temps suivant:kt+1 
       * 
       */
      /* L100: */
    }
  return 0;
}				/* icse1_ */
