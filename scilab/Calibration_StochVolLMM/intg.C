 #include "intg.h"

/*
based on SciLab's (Jean-Philippe Chancelier) and quadpack routines 
updated for Premia by Vadim Zherder, INRIA
*/ 
     const double epmach=1e-15, oflow=1e+100, uflow=1e-100;
      
      int iero;
      long int jupbnd;
     
      double epsr=1.0e-10;
      double epsa=1.0e-14;
      int limexp=50;
       
//-----------------------------------------------------------------------------------------------------//

void quarul(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double     *resasc);

void dqag0(double (*f)(double), double a, double b, double epsabs, double epsrel, double *result, double        *abserr, double *work, long int lwork, long int *iwork, long int liwork, int *ifail);

void dqags(double (*f)(double), double a, double b, double epsabs, double epsrel, double *alist, double         *blist,double *elist, double *rlist, long int limit, long int *iord, long int liord, double *result,        double *abserr, int *ier);

void order(long int limit, long int last, long int *maxerr, double *ermax, double *elist, long int *iord,        long int liord, long int nrmax);

void epsalg(long int n, double *epstab, double *result, double *abserr, double *res3la, long int *nres);

double min(double a, double b);

double max(double a, double b);

//****************************************************************************************************//
//------------------------------------INTG-----------------------------------------------------------//
//****************************************************************************************************//      

double * xconst, (*fmulti_origin)(double []);
int nconst;
void intg_multi(double a, double b, double (*f)(double []), double ea, double er, double *val, double *abserr, double x[], int n)
/* integral the function double f(double x[n]) with respect to x[0] in domain [a, b], while x[1]...x[n] is taken as given constants, modified by Xiao Wei, INRIA  */
{
 	nconst=n;
	xconst=(double *)calloc(nconst+1, sizeof (double));
	
	for (int i=1; i<=nconst; i++)
	{
		xconst[i]=x[i];
	}
	fmulti_origin=f;
	intg(a, b, f_new, ea, er, val, abserr);
	free(xconst);
}
extern double *xconst, (*fmulti_origin)(double []);
extern int nconst;
double f_new(double y)//from multifactors function to one factor function
{
	double *x=(double *)calloc(nconst+1, sizeof (double));
	x[0]=y;
	for (int i=1; i<=nconst; i++)
	{
		x[i]=xconst[i];
	}
	return double((*fmulti_origin)(x));
	free(x);
}


void intg(double a, double b, double (*f)(double), double ea, double er, double *val, double *abserr)
      {

      long int lw,liw;
      int ifail;
      double *work;
      long int *iwork;
      lw=3000;
      liw=3000/8+2;
        work=new double[lw+1];
	iwork=new long int[liw+1];
	
       dqag0(f,a,b,epsa,epsr,val,abserr, work, lw, iwork, liw, &ifail);
       
       delete[] work;
       delete[] iwork;

}
//****************************************************************************************************//
//--------------------------------QUARUL-----------------------------------------------------------//
//****************************************************************************************************//

void quarul(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double *resasc)
{
//
//     based on quadpack routine quarul
//     ******************************************************
//
//           purpose
//              to compute i = integral of f over (a,b), with error
//                             estimate
//                         j = integral of abs(f) over (a,b)
//
//           calling sequence
//              call quarul (f,a,b,result,abserr,resabs,resasc)
//
//           parameters
//              f      - function subprogram defining the integrand
//                       function f(x). the actual name for f needs
//                       to be declared e x t e r n a l in the
//                       calling program
//
//              a      - lower limit of integration
//
//              b      - upper limit of integration
//
//              result - approximation to the integral i.
//                       result is calculated by applying
//                       the 21-point gauss-kronrod rule
//                       (resk), obtained by optimal
//                       addition of abscissae to the
//                       10-point gauss  rule (resg).
//
//              abserr - estimate of the modulus of the
//                       absolute error, which should not
//                       exceed abs(i-result)
//              resabs - approximation to the integral j
//
//              resasc - approximation to the integral of
//                       abs(f-i/(b-a)) over (a,b)
//
//     ******************************************************
      double absc, centre, dhlgth, fc, fsum, fval1;
      double fval2,hlgth, resg, resk, reskh;
      long int j;
//     .. local arrays ..
      double fv1[11], fv2[11];// wg[11], wgk[11], xgk[11];
//            the abscissae and weights are given for the
//            interval (-1,1) . because of symmetry only the
//            positive abscissae and their corresponding
//            weights are given.
//            xgk    - abscissae of the 21-point gauss-kronrod rule
//                     xgk(2), xgk(4), .... abscissae of the 10-point
//                     gauss rule
//                     xgk(1), xgk(3), .... abscissae which
//                     are optimally added to the 10-point
//                     gauss rule
//            wgk    - weights of the 21-point gauss-kronrod rule
//            wg     - weights of the 10-point gauss rule,
//                     corresponding to the abscissae xgk(2),
//                     xgk(4), ... wg(1), wg(3), ... are set
//                     to zero.
//
     double xgk[12]={0.0, 0.99565716302580808073552728070,
     0.97390652851717172007796401210,
     0.93015749135570822600120718010,
     0.86506336668898451073209668840,
     0.78081772658641689706371757830,
     0.67940956829902440623432736510,
     0.56275713466860468333900009930,
     0.43339539412924719079926594320,
     0.29439286270146019813112660310,
     0.14887433898163121088482600110,0.0};
     double wgk[12]={0.0, 0.011694638867371874278064396060,
     0.032558162307964727478818972460,
     0.054755896574351996031381300240,
     0.075039674810919952767043140920,
     0.093125454583697605535065465080,
     0.10938715880229764189921059030,
     0.12349197626206585107795810980,
     0.13470921731147332592805400180,
     0.14277593857706008079709427310,
     0.14773910490133849137484151600,
     0.14944555400291690566493646840};
     double wg[11]={0.0, 0.0,0.066671344308688137593568809890,0.0,
     0.14945134915058059314577633970,0.0,
     0.21908636251598204399553493420,0.0,
     0.26926671930999635509122692160,0.0,
     0.29552422471475287017389299470};

//
//           list of major variables
//           ----------------------
//           centre - mid point of the interval
//           hlgth  - half length of the interval
//           absc   - abscissa
//           fval*  - function value
//           resg   - 10-point gauss formula
//           resk   - 21-point gauss-kronrod formula
//           reskh  - approximation to mean value of f over
//                    (a,b), i.e. to i/(b-a)
//
      centre = 0.5*(a+b);
      hlgth = 0.5*(b-a);
      dhlgth = fabs(hlgth);
//
//           compute the 21-point gauss-kronrod approximation to
//           the integral, and estimate the absolute error
//
      resg = 0.0;
      fc = f(centre);
      if(iero!=0) return;
      resk = wgk[11]*fc;
      *resabs = fabs(resk);

      for(j=1;j<=10;j++)
      {
         absc = hlgth*xgk[j];
         fval1 = f(centre-absc);
         if(iero!=0) {return;}
         fval2 = f(centre+absc);
         if(iero!=0) {return;}
         fv1[j] = fval1;
         fv2[j] = fval2;
         fsum = fval1 + fval2;
         resg = resg + wg[j]*fsum;
         resk = resk + wgk[j]*fsum;
         *resabs = *resabs + wgk[j]*(fabs(fval1)+fabs(fval2));
	 }

      reskh = resk*0.5;
      *resasc = wgk[11]*fabs(fc-reskh);

      for(j=1;j<=10;j++)
         {*resasc = *resasc + wgk[j]*(fabs(fv1[j]-reskh)+fabs(fv2[j]-reskh));}
      *result = resk*hlgth;
      *resabs = (*resabs)*dhlgth;
      *resasc = (*resasc)*dhlgth;
      *abserr = fabs((resk-resg)*hlgth);

      if( ((*resasc)!=0.0) && ((*abserr)!=0.0)) 
      	{
	*abserr =(*resasc)*min(1.0,pow((200.0*(*abserr)/(*resasc)),1.5));
	}
      if ((*resabs)>uflow/(50.0*epmach)) 
      {
	*abserr =max(epmach*(*resabs)*50.0,*abserr);
	}
      return;
}//quarul

//****************************************************************************************************//
//----------------------------DQAG-0-------------------------------------------------------------//
//****************************************************************************************************//

void dqag0(double (*f)(double), double a, double b, double epsabs, double epsrel, double *result, double *abserr, double *work, long int lwork, long int *iwork, long int liwork, int *ifail)
{
//
//     based on quadpack routine
//
//     dqag0 is a general purpose integrator which calculates
//     an approximation to the integral of a function over a finite
//     interval (a,b)
//
//     dqag0 itself is essentially a dummy routine whose function is to
//     partition the work arrays work and iwork for use by dqags.
//     work is partitioned into 4 arrays each of size lwork/4.
//     iwork is a single array in dqags.
//
//     .. scalar arguments ..
      
      int ibl, iel, ier, irl, limit;
//     .. subroutine references ..
//     dqags
//     ..
//     check that minimum workspace requirements are met
      if (lwork<4) 
      {
      	ier = 6;
      	*ifail = 1;
      	return;
      }
  
      if (liwork<lwork/8+2)
      {
      	ier = 6;
      	*ifail = 1;
      	return;
      }
//     limit = upper bound on number of subintervals
      limit = lwork/4;
//     set up base addresses for work arrays
      ibl = limit + 1;
      iel = limit + ibl;
      irl = limit + iel;
//     perform integration
      dqags(f, a, b, fabs(epsabs), fabs(epsrel), &(work[1]),&(work[ibl]), &(work[iel]), &(work[irl]), limit, iwork, liwork, result, abserr, &ier);
      
      if (ier!=0) {*ifail = 1;}
      else {*ifail = 0;}
      return;
 }//dqag0

//****************************************************************************************************//
//----------------------------------DQAGS----------------------------------------------
//****************************************************************************************************//
     void dqags(double (*f)(double), double a, double b, double epsabs, double epsrel, double *alist, double *blist,double *elist, double *rlist, long int limit, long int *iord, long int liord, double *result, double *abserr, int *ier)
     {
//
//     based on quadpack routine dqags (formerly qags)
//     **********************************************************
//
//        purpose
//           the routine calculates an approximation
//           /result/ to a given definite integral   i =
//           integral of /f/ over (a,b), hopefully
//           satisfying following claim for accuracy .
//           abs(i-result) .le. max(epsabs,epsrel*abs(i)).
//
//          calling sequence
//           call dqags (f,a,b,epsabs,epsrel,alist,blist,elist,
//                        rlist,limit,iord,liord,result,abserr,ier)
//
//        parameters
//            f      - function subprogram defining the integrand
//                     function f(x). the actual name for f
//                     needs to be declared e x t e r n a l
//                     in the driver program
//
//            a      - lower limit of integration
//
//            b      - upper limit of integration
//
//            epsabs - absolute accuracy requested
//
//            epsrel - relative accuracy requested
//
//            alist,blist,elist,rlist
//                   - work arrays (functions described below)
//
//            limit  - upper bound for number of subintervals
//
//            iord   - work array
//
//            liord  - length of iord (at least limit/2 + 2)
//
//            result - approximation to the integral
//
//            abserr - estimate of the modulus of the absolute error,
//                     which should equal or exceed abs(i-result)
//
//            ier    - ier   = 0 normal and reliable
//                             termination of the routine.
//                             it is assumed that the
//                             requested  accuracy has been
//                             achieved.
//                   - ier   .ne. 0 abnormal termination of
//                             the routine. the estimates
//                             for integral and error are
//                             less reliable. it is assumed
//                             that the  requested accuracy
//                             has not been achieved.
//                         = 1 maximum number of subdivisions allowed
//                             has been achieved. the user can
//                             allow more sub divisions by
//                             increasing the dimensions of the
//                             work arrays work and iwork.
//                             however, this may
//                             yield no  improvement, and it
//                             is rather advised to have a
//                             close look at the integrand,
//                             in order to determine the
//                             integration  difficulties. if
//                             the position of a local
//                             difficulty can be determined
//                             (i.e.  singularity,
//                             discontinuity within the
//                             interval) one will probably
//                             gain from  splitting up the
//                             interval at this point and
//                             calling the integrator on the
//                             sub-ranges. if possible, an
//                             appropriate special-purpose
//                             integrator should be used
//                             which is designed for
//                             handling the type  of
//                             difficulty involved.
//                         = 2 the occurrence of roundoff
//                             error is detected which
//                             prevents the requested
//                             tolerance  from being
//                             achieved. the error may be
//                             under-estimated.
//                         = 3 extremely bad integrand behaviour
//                             occurs at some interior points of the
//                             integration interval.
//                         = 4 it is presumed that the requested
//                             tolerance cannot be achieved,
//                             and that the returned result
//                             is the best which can be
//                             obtained.
//                         = 5 the integral is probably divergent, or
//                             slowly convergent. it must be noted
//                             that divergency can occur
//                             with any other value of ier.
//                         = -1 an error occurs during the evaluation of f
//     **********************************************************
    
//     .. local scalars ..
   double a1, a2, abseps, area12, area1, area2, area, b1;
     double b2,correc=0., defab1, defab2, defabs, dres, erlarg=0.,erlast;
     double errbnd, errmax, erro12, error1, error2, errsum,ertest=0.;
     double resabs, reseps, small=0.;
     long  int id, ierro, iroff1, iroff2, iroff3, k, ksgn, ktmin,last1;
     long int last=0, maxerr, nres, nrmax, numrl2;
     int extrap, noext;
     int flag, itsok=0, itsnotok=1;
//     .. local arrays ..
      double res3la[4], rlist2[53];
      
      iero=0;
//
//            list of major variables
//            -----------------------
//
//           alist     - list of left end-points of all subintervals
//                       considered up to now
//
//           blist     - list of right end-points of all subintervals
//                       considered up to now
//
//           rlist(i)  - approximation to the integral over
//                       (alist(i),blist(i))
//
//           rlist2    - array of dimension at least limexp+2
//                       containing the part of the epsilon table
//                       which is still needed for further
//                       computations
//
//           elist(i)  - error estimate applying to rlist(i)
//
//           maxerr    - pointer to the interval with largest error
//                       estimate
//
//           errmax    - elist(maxerr)
//
//           erlast    - error on the interval currently subdivided
//                       (before that subdivision has taken place)
//
//           area      - sum of the integrals over the subintervals
//
//           errsum    - sum of the errors over the subintervals
//
//           errbnd    - requested accuracy max(epsabs,epsrel*
//                       abs(result))
//
//           *****1    - variable for the left interval
//
//           *****2    - variable for the right interval
//
//           last      - index for subdivision
//
//           nres      - number of calls to the extrapolation routine
//
//           numrl2    - number of elements currently  in
//                       rlist2. if an appropriate
//                       approximation to the compounded
//                       integral has been obtained it is
//                       put in  rlist2(numrl2) after numrl2
//                       has been increased by one.
//
//           small     - length of the smallest interval considered
//                       up to now, multiplied by 1.5
//
//           erlarg    - sum of the errors over the intervals larger
//                       than the smallest interval
//                       considered up to now
//           extrap    - logical variable denoting that the
//                       routine is attempting to perform
//                       extrapolation.  i.e. before
//                       subdividing the smallest interval
//                       we try to decrease the value of
//                       erlarg
//           noext     - logical variable denoting that extrapolation
//                       is no longer allowed(/true/ value)
//
//           first approximation to the integral
//           -----------------------------------
//
      last1 = 1;
      *ier = 0;
      ierro = 0;
      quarul(f, a, b, result, abserr, &defabs, &resabs);
 
      if(iero>0){
         *ier=6;
         return;
	 }
	
//
//           test on accuracy
//
      dres = fabs(*result);
      errbnd = max(epsabs, epsrel*dres);
      if( (*abserr<=100*epmach*defabs)&&(*abserr>errbnd)){*ier = 2;}
      if( (limit<2) && (*abserr>errbnd) ) {*ier = 1;}
      if( (*ier!=0) || (*abserr<=errbnd) ) {
      	if (*ier>2) {*ier = *ier - 1;}
      	iord[1] = 4*last1;
      	return;
	}
//
//           initialization
//           --------------
//
      alist[1] = a;
      blist[1] = b;
      rlist[1] = *result;
      rlist2[1] = *result;
      errmax = *abserr;
      maxerr = 1;
      area = *result;
      errsum = *abserr;
      *abserr = oflow;
      nrmax = 1;
      nres = 0;
      numrl2 = 2;
      ktmin = 0;
      extrap = 0;
      noext = 0;
      iroff1 = 0;
      iroff2 = 0;
      iroff3 = 0;
      ksgn = -1;
      if (dres>=(1.0-50.0*epmach)*defabs) {ksgn = 1;}
//
//           main do-loop
//           ------------
//
      if (limit>=2)
      {
          for (last=2;last<=limit;last++)
	  {
//
//           bisect the subinterval with the nrmax-th largest
//           error estimate
//
               last1 = last;
	       a1 = alist[maxerr];
	       b1 = 0.5*(alist[maxerr]+blist[maxerr]);
	       a2 = b1;
	       b2 = blist[maxerr];
	       erlast = errmax;
	       quarul(f, a1, b1, &area1, &error1, &resabs, &defab1);
	       if(iero>0) {
	       	    *ier=6;
		    return;
		    }
	       quarul(f, a2, b2, &area2, &error2, &resabs, &defab2);
	       if(iero>0) {
	       	    *ier=6;
		    return;
		    }
//
//           improve previous approximation of integral
//           and error and test for accuracy
//
               area12 = area1 + area2;
	       erro12 = error1 + error2;
	       errsum = errsum + erro12 - errmax;
	       area = area + area12 - rlist[maxerr];
	       if( (defab1!=error1) && (defab2!=error2) ){
         	   if( (fabs(rlist[maxerr]-area12)<=0.00001*fabs(area12) )&& (erro12>=0.99*errmax))
		   {
         		if (extrap) {iroff2 = iroff2 + 1;}
         		else {iroff1 = iroff1 + 1;}
			}
          	    if ((last>10) || (erro12>errmax)) {iroff3 = iroff3 + 1;}
		    }
       
                rlist[maxerr] = area1;
		rlist[last] = area2;
		errbnd = max(epsabs,epsrel*fabs(area));
		if (errsum<=errbnd) {itsok=1; itsnotok=0;break;}
//           test for roundoff error and eventually
//           set error flag
//
                if ((iroff1+iroff2>=10) || (iroff3>=20)) {*ier = 2;}
		if (iroff2>=5) {ierro = 3;}
//
//           set error flag in the case that the number of interval
//            bisections exceeds /limit/
//
                if (last==limit) {*ier = 1;}
//
//           set error flag in the case of bad integrand behaviour
//           at interior points of integration range
//
                if (max(fabs(a1),fabs(b2))<=(1.1+100.0*epmach)*(fabs(a2)+1000.0*uflow))  {*ier = 4;}
		if (*ier!=0){break;}
//
//           append the newly-created intervals to the list
//
                if (error2>error1)
		{
			 alist[maxerr] = a2;
			 alist[last] = a1;
			 blist[last] = b1;
			 rlist[maxerr] = area2;
			 rlist[last] = area1;
			 elist[maxerr] = error2;
			 elist[last] = error1;
			 }
		else
		{
		         alist[last] = a2;
			 blist[maxerr] = b1;
			 blist[last] = b2;
			 elist[maxerr] = error1;
			 elist[last] = error2;
			 }
//
//           call subroutine order to maintain the
//           descending ordering in the list of error
//           estimates and select the subinterval with
//           nrmax-th largest error estimate (to be bisected
//           next)
//
                order(limit, last, &maxerr, &errmax, elist, iord,liord, nrmax);
		if (last==2) 
		{
		     small = fabs(b-a)*0.375;
		     erlarg = errsum;
		     ertest = errbnd;
		     rlist2[2] = area;
		     }
		if (noext) {continue;}
		erlarg = erlarg - erlast;
		if (fabs(b1-a1)>small) {erlarg = erlarg + erro12;}
		if (!extrap) 
		{
//
//           test whether the interval to be bisected next is the
//           smallest interval
//
                     if (fabs(blist[maxerr]-alist[maxerr])>small) {continue;}
		     extrap = 1;//.true.
		     nrmax = 2;
		     }
                if ((ierro!=3) && (erlarg>ertest)) 
		{
//
//           the smallest interval has the largest error.
//           before bisecting decrease the sum of the errors
//           over the larger intervals(erlarg) and perform
//           extrapolation
//
                     id = nrmax;
		     for (k=id;k<=jupbnd;k++)
		     {
		          maxerr = iord[nrmax];
			  errmax = elist[maxerr];
			  if (fabs(blist[maxerr]-alist[maxerr])>small) {break;}
			  nrmax = nrmax + 1;
			  }
	             if (fabs(blist[maxerr]-alist[maxerr])>small) {continue;}
		     }
//
//           perform extrapolation
//
                 numrl2 = numrl2 + 1;
		 rlist2[numrl2] = area;
		 epsalg(numrl2, rlist2, &reseps, &abseps, res3la, &nres);
		 ktmin = ktmin + 1;
		 if ((ktmin>5) && (*abserr<0.001*errsum)) {*ier = 5;}
		 if (abseps<*abserr) 
		 {
		     ktmin = 0;
		     *abserr = abseps;
		     *result = reseps;
		     correc = erlarg;
		     ertest = max(epsabs,epsrel*fabs(reseps));
		     if (*abserr<=ertest) {break;}
		     }
//
//           prepare  bisection of the smallest interval
//
                 if (numrl2==1) {noext = 1;}
		 if (*ier==5) {break;}
		 maxerr = iord[1];
		 errmax = elist[maxerr];
		 nrmax = 1;
		 extrap = 0;//.false.
		 small = small*0.5;
		 erlarg = errsum;
		 }//for last
//
//           set  final result and error estimate
//           ------------------------------------
//
        }//if limit>=2
 
        if(itsnotok)
	{
	    if (*abserr==oflow) {itsok=1;}
	    else 
	    {
	        flag=0;
             	if (*ier+ierro==0) {flag=1;}
        	else
		{
			if (ierro==3) {*abserr = *abserr + correc;}
        		if (*ier==0) {*ier = 3;}
        		if( (*result!=0.0) && (area!=0.0)) 
                        {
		        	if (*abserr/fabs(*result)>errsum/fabs(area)) {itsok=1;flag=0;}
				else {flag=1;}
				}
			else
			{
				if (*abserr>errsum) {itsok=1; flag=0;}
				else
        			{
                                     if (area==0.0) {itsok=0;flag=0;}
				     else {flag=1;}
				     }
				}
			}
		
//           test on divergency
//
	     if(flag)
	     {
	        if ((ksgn==-1) && (max(fabs(*result),fabs(area))<=defabs*0.01)) {itsok=0;}
		else 
		{
			if( (0.01>(*result/area) || (*result/area)>100.0) || (errsum>fabs(area)) ) {*ier = 6;}
			}
		itsok=0;
		}
	   }
	} // if itsnotok
//
//           compute global integral sum
//
 // 280 
      if(itsok)
      {
             *result = 0.0;
	     for (k=1;k<=last;k++) {*result = *result + rlist[k];}
	     *abserr = errsum;
      }

      if (*ier>2) {*ier = *ier - 1;}
      iord[1] = 4*last1;
      return;
 }//end dqags
 
//****************************************************************************************************// 
//------------------------------ORDER -------------------------------------------------------------//
//****************************************************************************************************//

void order(long int limit, long int last, long int *maxerr, double *ermax, double *elist, long int *iord, long int liord, long int nrmax)
{
//
//     based on quadpack routine order
//     ******************************************************
//
//           purpose
//              this routine maintains the descending ordering
//              in the list of the local error estimates
//              resulting from the interval subdivision
//              process. at each call two error estimates
//              are inserted using the sequential search
//              method . top-down for the largest error
//              estimate,  bottom-up for the smallest error
//              estimate.
//
//           calling sequence
//              call order
//              (limit,last,maxerr,ermax,elist,iord,liord,nrmax)
//
//             parameters (meaning at output)
//              limit  - maximum number of error estimates the list
//                       can contain
//
//              last   - number of error estimates currently
//                       in the list. elist(last) contains
//                       the smallest error estimate.
//
//              maxerr - maxerr points to the nrmax-th largest error
//                       estimate currently in the list.
//
//              ermax  - nrmax-th largest error estimate
//                       ermax = elist(maxerr)
//
//              elist  - array of dimension last containing
//                       the error estimates
//
//              iord   - array containing pointers to elist so
//                       that iord(1) points to the largest
//                       error estimate,...,iord(last) to the
//                       smallest error estimate
//
//              liord  - dimension of iord
//
//              nrmax  - maxerr = iord(nrmax)
//
//     ******************************************************
//
//     .. local scalars ..
      double errmax, errmin;
      long int i, ibeg, ido, isucc, j, jbnd, k;
      int flag;
//            check whether the list contains more than
//            two error estimates
//
      if (last<=2)
      { 
      	iord[1] = 1;
      	iord[2] = 2;
	}
      else
      {
//
//           this part of the routine is only executed
//           if, due to a difficult integrand, subdivision
//           increased the error estimate. in the normal case
//           the insert procedure should start after the
//           nrmax-th largest error estimate.
//
             errmax = elist[*maxerr];
	     if (nrmax!=1)
	     {
	          ido = nrmax - 1;
		  for(i=1;i<=ido;i++)
		  {
		      isucc = iord[nrmax-1];
		      if (errmax<=elist[isucc]) {break;}
		      iord[nrmax] = isucc;
		      nrmax = nrmax - 1;
		      }
		  }
//
//           compute the number of elements in the list to
//           be maintained in descending order. this number
//           depends on the number of subdivisions still
//           allowed
             flag=0;
	     jupbnd = last;
	     if (last>=(limit/2+2)) {jupbnd = limit + 3 - last;}
	     errmin = elist[last];
//
//           insert errmax by traversing the list top-down
//           starting comparison from the element
//           elist(iord(nrmax+1))
//
             jbnd = jupbnd - 1;
	     ibeg = nrmax + 1;
	     if (ibeg<=jbnd)
	     {
	         for(i=ibeg;i<=jbnd;i++)
		 {
		     isucc = iord[i];
		     if (errmax>=elist[isucc])
		     {
		     	 flag=1;
			 iord[i-1] = *maxerr;
			 k = jbnd;
			 for(j=i;j<=jbnd;j++)
			 {
			      isucc = iord[k];
			      if (errmin<elist[isucc]) {iord[k+1] = last;break;}
			      iord[k+1] = isucc;
			      k = k - 1;
			      }
			 if (errmin>=elist[isucc]) {iord[i] = last;}
			 break;
			 }
		     iord[i-1] = isucc;
		     }//for i
		 }//if ibeg

             if(!flag)
	     {
	         iord[jbnd] = *maxerr;
		 iord[jupbnd] = last;
		 }
        }//else
	
      *maxerr = iord[nrmax];
      *ermax = elist[*maxerr];
      return;
}// order

//****************************************************************************************************//
//-------------------------EPSALG----------------------------------------------------------//
//****************************************************************************************************//

void epsalg(long int n, double *epstab, double *result, double *abserr, double *res3la, long int *nres)
{
//
//     based on quadpack routine epsalg
//     ******************************************************
//
//           purpose
//              the routine transforms a given sequence of
//              approximations, by means of the epsilon
//              algorithm of p. wynn.
//
//              an estimate of the absolute error is also given.
//              the condensed epsilon table is computed. only those
//              elements needed for the computation of the
//              next diagonal are preserved.
//
//           calling sequence
//              call epsalg (n,epstab,result,abserr,res3la,nres)
//
//           parameters
//              n      - epstab(n) contains the new element in the
//                       first column of the epsilon table.
//
//              epstab - one dimensional array containing the
//                       elements of the two lower diagonals of
//                       the triangular epsilon table.
//                       the elements are numbered starting at the
//                       right-hand corner of the triangle.
//                       the dimension should be at least n+2.
//
//              result - resulting approximation to the integral
//
//              abserr - estimate of the absolute error computed from
//                       result and the 3 previous /results/
//
//              res3la - array containing the last 3 /results/
//
//              nres   - number of calls to the routine
//                       (should be zero at first call)
//
//     ******************************************************
      
//     .. local scalars ..
      double delta1, delta2, delta3, e0, e1, e1abs, e2, e3, epsinf;
      double err1, err2, err3, error, res, ss, tol1, tol2, tol3;
      long int i, ib2, ib, ie, ind, k1, k2, k3, newelm, num;
//             -------------------------
//            /limexp/ is the maximum number of elements the epsilon
//            table can contain. if this number is reached, the upper
//            diagonal of the epsilon table is deleted.
//
      
//           list of major variables
//           -----------------------
//           e0     - the 4 elements on which the
//           e1       computation of a new element in
//           e2       the epsilon table is based
//           e3                 e0
//                        e3    e1    new
//                              e2
//           newelm - number of elements to be computed in the new
//                    diagonal
//           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
//           result - the element in the new diagonal with least
//                    error
//
      *nres = *nres + 1;
      *abserr = oflow;
      *result = epstab[n];
      if (n<3) 
      {
        	*abserr = max(*abserr,5.0*epmach*fabs(*result));
	        return;
		}
      epstab[n+2] = epstab[n];
      newelm = (n-1)/2;
      epstab[n] = oflow;
      num = n;
      k1 = n;

      for(i=1;i<=newelm;i++)
      {
         k2 = k1 - 1;
         k3 = k1 - 2;
         res = epstab[k1+2];
         e0 = epstab[k3];
         e1 = epstab[k2];
         e2 = res;
         e1abs = fabs(e1);
         delta2 = e2 - e1;
         err2 = fabs(delta2);
         tol2 = max(fabs(e2),e1abs)*epmach;
         delta3 = e1 - e0;
         err3 = fabs(delta3);
         tol3 = max(e1abs,fabs(e0))*epmach;
         if ((err2<=tol2) && (err3<=tol3)) 
	 {
//
//           if e0, e1 and e2 are equal to within machine
//           accuracy, convergence is assumed
//           result = e2
//           abserr = abs(e1-e0)+abs(e2-e1)
//
         	 *result = res;
		 *abserr = err2 + err3;
		 *abserr = max(*abserr,5.0*epmach*fabs(*result));
		 return;
		 }
     	 e3 = epstab[k1];
         epstab[k1] = e1;
         delta1 = e1 - e3;
         err1 = fabs(delta1);
         tol1 = max(e1abs,fabs(e3))*epmach;
//
//           if two elements are very close to each other, omit
//           a part of the table by adjusting the value of n
//
         if ((err1<tol1) || (err2<tol2) || (err3<tol3)) 
	 {
	       n = i + i - 1;
	       break;
	       }
         ss = 1.0/delta1 + 1.0/delta2 - 1.0/delta3;
         epsinf = fabs(ss*e1);
//
//           test to detect irregular behaviour in the table, and
//           eventually omit a part of the table adjusting the value
//           of n
//
         if (epsinf<=0.0001)
	 { 
	       n = i + i - 1;
	       break;
	       }
//
//           compute a new element and eventually adjust
//           the value of result
//
         res = e1 + 1.0/ss;
         epstab[k1] = res;
         k1 = k1 - 2;
         error = err2 + fabs(res-e2) + err3;
         if (error>=(*abserr)) {continue;}
         *abserr = error;
         *result = res;
      }//for
//           shift the table
//
      if (n==limexp) {n = 2*(limexp/2) - 1;}
      ib = 1;
      if ((num/2)*2==num) {ib = 2;}
      ie = newelm + 1;

      for(i=1;i<=ie;i++)
      {
         ib2 = ib + 2;
         epstab[ib] = epstab[ib2];
         ib = ib2;
	 }

      if (num!=n) 
      {
         ind = num - n + 1;
	 for(i=1;i<=n;i++)
	 {
	    epstab[i] = epstab[ind];
	    ind = ind + 1;
	    }
	 }

      if (*nres<4)
      {
         res3la[*nres] = *result;
	 *abserr = oflow;
	 }
      else
      {
//           compute error estimate
//
          *abserr = fabs(*result-res3la[3]) + fabs(*result-res3la[2]) +fabs(*result-res3la[1]);
	  res3la[1] = res3la[2];
	  res3la[2] = res3la[3];
	  res3la[3] = *result;
	  }
	  
      *abserr = max(*abserr,5.0*epmach*fabs(*result));
      return;
 }//epsalg
//****************************************************************************************************//
double min(double a, double b)
{
	double res=a>b?b:a;
	return res;
} 
 double max(double a, double b)
 {
 	double res=a>b?a:b;
	return res;
 }

