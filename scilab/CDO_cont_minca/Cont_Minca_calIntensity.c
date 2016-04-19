#include"time.h"
#include"pnl_vector.h"
#include"pnl_random.h"
#include"pnl_specfun.h"
#include"pnl_optim.h"
#include"pnl_matrix.h"
#include"pnl_mathtools.h"

//%%%%%%%%%%%%%%%%%%%%%%%% Optimal u(t,L) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void nl_constraints_func(const PnlVect *x,PnlVect *fx,void *params)
{
  pnl_vect_resize(fx,1);
  pnl_vect_set_double(fx,1.0);
}
//This function returns the optimal value of -V(0,0,mu0) via dynamic programming.
double optu(const PnlVect *x,void *param)
{
	double fx,recov,dz,dz1,Bdz,sum,sum1,sum2,*params;
	int i,j,numName,numTran,numMat,n,iMat,k;
	PnlVect *temp_vect1,*temp_vect2,*attpts,*numPeriod,*B,*U;
  	PnlMat *mu0,*muu,*mud,*mu,*Phicoeff,*Phicoeff1,*PhicoeffJ,*Phi,*temp_mat,*ePhi,*u
		,*ub,*lb,*s,*c,*Vcoeff_dz,*Vcoeff_dz1;

        params=(double *)param;
	attpts=pnl_vect_create(pnl_iround(params[0]));
	ub=pnl_mat_create(pnl_iround(params[1]),pnl_iround(params[2]));
	lb=pnl_mat_create(pnl_iround(params[1]),pnl_iround(params[2]));
	numPeriod=pnl_vect_create(pnl_iround(params[9]));
	B=pnl_vect_create(pnl_iround(params[3]));
	U=pnl_vect_create(pnl_iround(params[4]));
	s=pnl_mat_create(pnl_iround(params[7]),pnl_iround(params[8]));
	c=pnl_mat_create(pnl_iround(params[0]),pnl_iround(params[5]));
	Vcoeff_dz=pnl_mat_create(pnl_iround(params[6]),pnl_iround(params[6]));
	Vcoeff_dz1=pnl_mat_create(pnl_iround(params[6]),pnl_iround(params[6]));

    n=14;
    numMat=pnl_iround(params[9]);
	numName=c->n-1;
    numTran=attpts->size;
    recov=params[10];
	dz=params[11];
    dz1=params[12];
	Bdz=params[13];

	for(i=0;i<attpts->size;i++)
    {
	   pnl_vect_set(attpts,i,params[i+n]);
	}
	for(i=0;i<numPeriod->size;i++)
    {
	   pnl_vect_set(numPeriod,i,params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+i]);
	}
	for(i=0;i<B->size;i++)
    {
           pnl_vect_set(B,i,params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+numPeriod->size+i]);
	}
	for(i=0;i<ub->m;i++)
        for(j=0;j<ub->n;j++)
	   {
	      pnl_mat_set(lb,i,j,params[n+attpts->size+ub->n*i+j]);
	      pnl_mat_set(ub,i,j,params[n+attpts->size+ub->n*ub->m+ub->n*i+j]);
	   }
	for(i=0;i<s->m;i++)
        for(j=0;j<s->n;j++)
	   {
	      pnl_mat_set(s,i,j,params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+numPeriod->size+B->size+s->n*i+j]);
	   }
	for(i=0;i<c->m;i++)
        for(j=0;j<c->n;j++)
	   {
	      pnl_mat_set(c,i,j,params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+numPeriod->size+B->size+s->n*s->m+c->n*i+j]);
	   }
	for(i=0;i<Vcoeff_dz->m;i++)
        for(j=0;j<Vcoeff_dz->n;j++)
	   {
	      pnl_mat_set(Vcoeff_dz,i,j,params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+numPeriod->size+B->size+s->n*s->m+c->n*c->m+Vcoeff_dz->n*i+j]);
	      pnl_mat_set(Vcoeff_dz1,i,j,params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+numPeriod->size+B->size+s->n*s->m+c->n*c->m+Vcoeff_dz->n*Vcoeff_dz->m+i*Vcoeff_dz->n+j]);
	   }
	for(i=0;i<U->size;i++)
    {
	   pnl_vect_set(U,i,params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+numPeriod->size+B->size+s->n*s->m+c->n*c->m+Vcoeff_dz->n*Vcoeff_dz->m+Vcoeff_dz->n*Vcoeff_dz->m+i]);
	}
    mu0=pnl_mat_create(pnl_iround(params[7]),pnl_iround(params[8]));
    for(j=0;j<mu0->n;j++)
  	   pnl_mat_set(mu0,0,j,pnl_vect_get(x,j));
    for(j=0;j<mu0->n;j++)
  	   pnl_mat_set(mu0,1,j,pnl_vect_get(x,mu0->n+j)/pnl_vect_get(attpts,0));
	for(i=2;i<mu0->m;i++)
           for(j=0;j<mu0->n;j++)
  				pnl_mat_set(mu0,i,j,pnl_vect_get(x,i*mu0->n+j)/(pnl_vect_get(attpts,i-1)-pnl_vect_get(attpts,i-2)));
	sum=0;
     for(i=0;i<ub->m;i++)
	   for(j=0;j<ub->n;j++)
		sum+=fabs(pnl_mat_get(ub,i,j))+fabs(pnl_mat_get(lb,i,j));
	muu=pnl_mat_create_from_double(mu0->m,mu0->n,0.0);
	mud=pnl_mat_create_from_double(mu0->m,mu0->n,0.0);
	mu=pnl_mat_create_from_double(mu0->m,mu0->n,0.0);
	if (sum> 0)
	{
       for(i=0;i<mu0->m;i++)
	     for(j=0;j<mu0->n;j++)
	     {
			pnl_mat_set(mu0,i,j,pnl_mat_get(mu0,i,j)*(pnl_mat_get(s,i,j)>0)); //set mu0=0 for s<0
	     }
       for(i=0;i<muu->m;i++)
	     for(j=0;j<muu->n && j<numMat;j++)
	     {
			pnl_mat_set(muu,i,j,pnl_mat_get(mu0,i,j));
	     }
       for(i=0;i<mud->m;i++)
	     for(j=0;j<mud->n && j+numMat<mud->n;j++)
	     {
			pnl_mat_set(mud,i,j,pnl_mat_get(mu0,i,j+numMat));
	     }
	   pnl_mat_clone(mu,mud);
       pnl_mat_minus_mat(mu,muu);
	}
	else
	{
	  pnl_mat_clone(muu,mu0);
	  pnl_mat_clone(mud,mu0);
      for(i=0;i<mu->m;i++)
	     for(j=0;j<mu->n;j++)
	     {
	    	pnl_mat_set(mu,i,j,pnl_mat_get(mu0,i,j)*(pnl_mat_get(s,i,j)>0)); //set mu=0 for s<0
	     }
	}

	Phicoeff1 = pnl_mat_create_from_double(numTran,numMat,0.0);                                                  //coefficients of Phi at t1
	Phicoeff = pnl_mat_create_from_double(numTran,numMat,0.0);                                                   //coefficients of Phi at j
	PhicoeffJ = pnl_mat_create_from_double(numTran,numMat,0.0);                                                  //coefficients of Phi at T

	for(i=0;i<numTran-1;i++)
          for(j=0;j<numMat;j++)
             pnl_mat_set(Phicoeff1,i,j,pnl_mat_get(mu,i+1,j)*(1+pnl_mat_get(s,i+1,j)*dz1-Bdz)-pnl_mat_get(mu,i+2,j)*(1+pnl_mat_get(s,i+2,j)*dz1-Bdz));
	for(i=0;i<numTran-1;i++)
          for(j=0;j<numMat;j++)
             pnl_mat_set(Phicoeff,i,j,pnl_mat_get(mu,i+1,j)*(1+pnl_mat_get(s,i+1,j)*dz-Bdz)-pnl_mat_get(mu,i+2,j)*(1+pnl_mat_get(s,i+2,j)*dz-Bdz));
	for(i=0;i<numTran-1;i++)
          for(j=0;j<numMat;j++)
             pnl_mat_set(PhicoeffJ,i,j,pnl_mat_get(mu,i+1,j)*(1+pnl_mat_get(s,i+1,j)*dz)-pnl_mat_get(mu,i+2,j)*(1+pnl_mat_get(s,i+2,j)*dz));

    for(j=0;j<numMat;j++)
	   pnl_mat_set(Phicoeff1,numTran-1,j,pnl_mat_get(mu,0,j)*(pnl_mat_get(s,0,j)*dz1/(1-recov)+1-Bdz)+pnl_mat_get(mu,mu->m-1,j)*(1+pnl_mat_get(s,s->m-1,j)*dz1-Bdz));
    for(j=0;j<numMat;j++)
       pnl_mat_set(Phicoeff,numTran-1,j,pnl_mat_get(mu,0,j)*(pnl_mat_get(s,0,j)*dz/(1-recov)+1-Bdz)+pnl_mat_get(mu,mu->m-1,j)*(1+pnl_mat_get(s,s->m-1,j)*dz-Bdz));
    for(j=0;j<numMat;j++)
       pnl_mat_set(PhicoeffJ,numTran-1,j,pnl_mat_get(mu,0,j)*(pnl_mat_get(s,0,j)*dz/(1-recov)+1)+pnl_mat_get(mu,mu->m-1,j)*(1+pnl_mat_get(s,s->m-1,j)*dz));
	Phi=pnl_mat_create_from_double(pnl_iround(pnl_vect_get(numPeriod,numPeriod->size-1)+1)*numMat,(numName+1),0.0); //compute (negative) Phi for all maturities

	temp_vect1=pnl_vect_create(PhicoeffJ->m);
	temp_vect2=pnl_vect_create(numName+1);
	temp_mat=pnl_mat_create(c->n,c->m);
	pnl_mat_tr(temp_mat,c);
	iMat = numMat-1;
	while (iMat >= 0)
	{
		pnl_mat_get_col(temp_vect1,PhicoeffJ,iMat);
		temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);
		for(j=0;j<Phi->n;j++)
			pnl_mat_set(Phi,iMat*(Phi->m/numMat)+pnl_iround(pnl_vect_get(numPeriod,iMat)),j,-pnl_vect_get(B,pnl_iround(pnl_vect_get(numPeriod,iMat))-1)
			*pnl_vect_get(temp_vect2,j));
		pnl_mat_get_col(temp_vect1,Phicoeff1,iMat);
		temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);
		for (i =pnl_iround(pnl_vect_get(numPeriod,iMat))-1;i>1;i--)
			for(j=0;j<Phi->n;j++)
				pnl_mat_set(Phi,iMat*(Phi->m/numMat)+i,j,-pnl_vect_get(B,i-1)*pnl_vect_get(temp_vect2,j));
		pnl_mat_get_col(temp_vect1,Phicoeff1,iMat);
		temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);
		for(j=0;j<Phi->n;j++)
			pnl_mat_set(Phi,iMat*(Phi->m/numMat)+1,j,-pnl_vect_get(B,0)*pnl_vect_get(temp_vect2,j));

		sum1=0;
		for(i=0;i<mu->m-2;i++)
			sum1+=pnl_mat_get(mu,i+2,iMat)*(pnl_vect_get(attpts,i+1)-pnl_vect_get(attpts,i));
		sum2=0;
		for(i=1;i<pnl_iround(pnl_vect_get(numPeriod,iMat));i++)
			sum2+=pnl_vect_get(B,i);

		for(j=0;j<Phi->n;j++)
			pnl_mat_set(Phi,iMat*(Phi->m/numMat),j,-(-pnl_vect_get(B,0)*(sum1+pnl_mat_get(mu,1,iMat)*pnl_vect_get(attpts,0)+pnl_mat_get(mu,0,iMat))
			+pnl_mat_get(mu,1,iMat)*pnl_vect_get(U,iMat)*pnl_vect_get(attpts,0)-pnl_mat_get(mu,0,iMat)*pnl_mat_get(s,0,iMat)*recov/(1-recov)
			*(pnl_vect_get(B,0)*dz1+sum2)*dz));

		iMat = iMat - 1;
	}
	sum=0;
    for(i=0;i<ub->m;i++)
	   for(j=0;j<ub->n;j++)
		sum+=pnl_mat_get(ub,i,j)*pnl_mat_get(muu,i,j)-pnl_mat_get(lb,i,j)*pnl_mat_get(mud,i,j);
	for(i=1;i<=numMat;i++)
		for(j=0;j<Phi->n;j++)
			pnl_mat_set(Phi,i*(Phi->m/numMat)-1,j,pnl_mat_get(Phi,i*(Phi->m/numMat)-1,j)+pnl_mat_get(Phi,(i-1)*(Phi->m/numMat),j)+sum);
	ePhi=pnl_mat_create(Phi->m/numMat,numName+1);
    for(i=0;i<ePhi->m;i++)
		for(j=0;j<ePhi->n;j++)
		{
			sum=0;
			for(k=0;k<numMat;k++)
				sum+=pnl_mat_get(Phi,k*(Phi->m/numMat)+i,j);
			pnl_mat_set(ePhi,i,j,exp(sum));//exp of sum of Phi in different maturities
		}
	u=pnl_mat_create(pnl_iround(pnl_vect_get(numPeriod,numPeriod->size-1))+1,ePhi->n);
	//dynamic programming
	pnl_vect_resize(temp_vect1,u->n);
	pnl_vect_resize(temp_vect2,u->n);
	pnl_mat_resize(temp_mat,Vcoeff_dz->n,Vcoeff_dz->m);
	pnl_mat_tr(temp_mat,Vcoeff_dz);

	pnl_mat_get_row(temp_vect1,ePhi,ePhi->m-1);
	pnl_mat_set_row(u,temp_vect1,u->m-1);

	temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);
	pnl_mat_set_row(u,temp_vect2,u->m-2);
	for (i=u->m-3;i>0;i--)
	{
		for(j=0;j<u->n;j++)
			pnl_vect_set(temp_vect1,j,(pnl_mat_get(u,i+1,j)*pnl_mat_get(ePhi,i+1,j)));
		temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);
		pnl_mat_set_row(u,temp_vect2,i);
	}
	pnl_mat_tr(temp_mat,Vcoeff_dz1);
	for(j=0;j<u->n;j++)
		pnl_vect_set(temp_vect1,j,(pnl_mat_get(u,1,j)*pnl_mat_get(ePhi,1,j)));

	temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);
	pnl_mat_set_row(u,temp_vect2,0);

	fx= log(pnl_mat_get(u,0,0));
	if (!pnl_isfinite(fx))//invalid function may affect calibration results
	{
		//printf("Warning: Invalid objective function. Calibration results maybe inaccurate.\n");
		fx = 1e6;
	}

	//FREE
	pnl_vect_free(&attpts);
	pnl_mat_free(&ub);
	pnl_mat_free(&lb);
	pnl_mat_free(&Vcoeff_dz);
	pnl_mat_free(&Vcoeff_dz1);
	pnl_mat_free(&s);
	pnl_mat_free(&c);
	pnl_vect_free(&B);
	pnl_vect_free(&U);

	pnl_mat_free(&mu0);
	pnl_mat_free(&mu);
	pnl_mat_free(&mud);
	pnl_mat_free(&muu);
	pnl_mat_free(&Phicoeff1);
	pnl_mat_free(&Phicoeff);
	pnl_mat_free(&PhicoeffJ);
	pnl_mat_free(&Phi);
	pnl_mat_free(&temp_mat);
	pnl_mat_free(&ePhi);
	pnl_mat_free(&u);
	pnl_vect_free(&temp_vect1);
	pnl_vect_free(&temp_vect2);
    pnl_vect_free(&numPeriod);

	return fx;
}
/*
% This function is to recover default intensities implied by CDO quotes
% via intensity control. (see Cont and Minca (2008) for details)
%
% Remarks:
% - Inputs accept index and tranche spreads for multiple maturities
% - Flexible time-to-maturity (not necessarily multiple of 0.25)
% - Arbitrary pior distribution (time-homogeneous intensity function)
% - Inequality constraints on mark-to-market value
% - Allow to ignore calibration results of certain spreads by setting the
% spread inputs to be negative
%
% Inputs
% T: time to maturity vector
% mktsprd: market spreads [index sprd(bps),equity upfront(%),2nd tranche sprd(bps),...]'
% attpts: attachment points
% lam0: dependence of prior intensity function on number of defaults
% ub: upper bounds on mark-to-market values
% lb: lower bounds on mark-to-market values
% mu0: initial value of mu for optimization
% muub: a common absolute bounds on mu (>1e6 for unbounded optimization)
%
% Outputs
% lambda: calibrated default intensity surface
*/
int calIntensity(double r,int numName,double dz,double dt,double recov,PnlVect *T0,PnlMat *mktsprd,PnlVect *attpts,PnlVect *lam0,PnlMat *ub,PnlMat *lb, PnlMat *mu00,double muub,double tolerance,int iter_max, PnlMat* lambda)
{
	double dz1,LGD,Bdz,sum,*optu_params,sum1,sum2;
	int i,j,numTran,numMat,print_inner_steps,result_optu,n,iMat,k;
	PnlVect *numPeriod,*B,*U,*lower_bounds,*upper_bounds,*mu_vect,*mu_vect0,*temp_vect1,*temp_vect2,*T;
	PnlMat *c,*Vcoeff_dz,*Vcoeff_dz1,*Vcoeff_dt,*temp1,*temp2,*s,*muu,*mud,*mu,*mu0,
		   *Phicoeff,*Phicoeff1,*PhicoeffJ,*Phi,*temp_mat,*ePhi,*u;
	PnlRnFuncR *optu_func;
   	PnlRnFuncRm *optu_func_grad,*nl_constraints;

	mu0=pnl_mat_create(mu00->m,mu00->n);
	pnl_mat_clone(mu0,mu00);
	mu_vect=pnl_vect_create(mu0->m*mu0->n);
	mu_vect0=pnl_vect_create(mu0->m*mu0->n);
    for(i=0;i<mu_vect->size;i++)
	   pnl_vect_set(mu_vect0,i,pnl_mat_get(mu0,(i/mu0->n),i%mu0->n));
	numPeriod=pnl_vect_create(T0->size);
	T=pnl_vect_create(T0->size);

	optu_func=(PnlRnFuncR *)malloc(sizeof(PnlRnFuncR));
    nl_constraints=(PnlRnFuncRm *)malloc(sizeof(PnlRnFuncRm));
	//parameters --------------------------------------------------------------
	LGD = (1-recov)/numName;//LGD
	for(i=0;i<T->size;i++)
	  pnl_vect_set(numPeriod,i,pnl_iround(trunc(pnl_vect_get(T0,i)/dz)));//no. of payment times
	dz1 = pnl_vect_get(T0,0)-trunc(pnl_vect_get(T0,0)/dz)*dz;//time-to-1st-payment
	if(dz1 == 0)
		dz1 = dz;
	else
		dz1 = floor(dz1/dt) * dt;//round time-to-1st-payment in scale of dt

	pnl_vect_clone(T,numPeriod);
	pnl_vect_minus_double(T,1);
	pnl_vect_mult_double(T,dz);
	pnl_vect_plus_double(T,dz1);
	numTran = attpts->size;                                                             //number of tranches
	numMat = T->size;                                                                   //number of maturities
	//settings for calibration ------------------------------------------------
	Bdz = exp(-r*dz);//1-period discount factor with time interval dz
	B=pnl_vect_create(pnl_iround(((pnl_vect_get(T,T->size-1)-dz1)/dz))+1);
	for(i=0;i<B->size;i++)
	  pnl_vect_set(B,i,exp(-r *(dz1+i*dz)));//discount factors
	c=pnl_mat_create(attpts->size,numName+1);
	for(i=0;i<c->m;i++)
		for(j=0;j<c->n;j++)
			pnl_mat_set(c,i,j,MAX(0,pnl_vect_get(attpts,i)-LGD*j));//E(K-L)^+
	s=pnl_mat_create(mktsprd->m,mktsprd->n);
	pnl_mat_clone(s,mktsprd);
	pnl_mat_div_double(s,10000.0);
	for(j=0;j<s->n;j++)
		pnl_mat_set(s,1,j,0.05);//s: sprds
	U=pnl_vect_create(mktsprd->n);//U: equity upfront
	pnl_mat_get_row(U,mktsprd,1);
	pnl_vect_div_double(U,100.0);

	Vcoeff_dz=pnl_mat_create_from_double(lam0->size,lam0->size,0.0);
	Vcoeff_dz1=pnl_mat_create_from_double(lam0->size,lam0->size,0.0);
	Vcoeff_dt=pnl_mat_create_from_double(lam0->size,lam0->size,0.0);
	temp1=pnl_mat_create_from_double(lam0->size,lam0->size,0.0);
	temp2=pnl_mat_create(lam0->size,lam0->size);
	for(i=0;i<lam0->size-1;i++)
	{
		pnl_mat_set(temp1,i,i,-pnl_vect_get(lam0,i));
		pnl_mat_set(temp1,i,i+1,pnl_vect_get(lam0,i));
	}
	pnl_mat_set(temp1,lam0->size-1,lam0->size-1,-pnl_vect_get(lam0,i));
	pnl_mat_sq_transpose(temp1);
	pnl_mat_mult_double(temp1,dz);
	pnl_mat_exp(temp2,temp1);
	for(i=0;i<Vcoeff_dz->m;i++)
		for(j=0;j<Vcoeff_dz->n;j++)
			pnl_mat_set(Vcoeff_dz,i,j,MAX(0,pnl_mat_get(temp2,i,j)));
	pnl_mat_mult_double(temp1,dz1/dz);
	pnl_mat_exp(temp2,temp1);
	for(i=0;i<Vcoeff_dz->m;i++)
		for(j=0;j<Vcoeff_dz->n;j++)
			pnl_mat_set(Vcoeff_dz1,i,j,MAX(0,pnl_mat_get(temp2,i,j)));
	pnl_mat_mult_double(temp1,dt/dz1);
	pnl_mat_exp(temp2,temp1);
	for(i=0;i<Vcoeff_dz->m;i++)
		for(j=0;j<Vcoeff_dz->n;j++)
			pnl_mat_set(Vcoeff_dt,i,j,MAX(0,pnl_mat_get(temp2,i,j)));
	//calibration -------------------------------------------------------------
    optu_func->function=optu;
	optu_func_grad=NULL;
    nl_constraints->function=nl_constraints_func;
	nl_constraints->params=NULL;
    n=14;
    optu_params=(double *)malloc((n+attpts->size+ub->n*ub->m+ub->n*ub->m+numPeriod->size+B->size+s->n*s->m+c->n*c->m
		+Vcoeff_dz->n*Vcoeff_dz->m+Vcoeff_dz->n*Vcoeff_dz->m+U->size)*sizeof(double));
    optu_params[0]=attpts->size;
	optu_params[1]=ub->m;
	optu_params[2]=ub->n;
	optu_params[3]=B->size;
	optu_params[4]=U->size;
	optu_params[5]=c->n;
	optu_params[6]=lam0->size;
	optu_params[7]=mu0->m;
	optu_params[8]=mu0->n;
    optu_params[9]=T->size;
    optu_params[10]=recov;
	optu_params[11]=dz;
    optu_params[12]=dz1;
	optu_params[13]=Bdz;
	for(i=0;i<attpts->size;i++)
	   optu_params[i+n]=pnl_vect_get(attpts,i);
	for(i=0;i<ub->m;i++)
	  for(j=0;j<ub->n;j++)
	  {
		optu_params[n+attpts->size+ub->m*i+j]=pnl_mat_get(lb,i,j);
		optu_params[n+attpts->size+ub->m*ub->n+ub->m*i+j]=pnl_mat_get(ub,i,j);
	  }
	for(i=0;i<numPeriod->size;i++)
    {
	   optu_params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+i]=pnl_vect_get(numPeriod,i);
	}
	for(i=0;i<B->size;i++)
    {
        optu_params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+numPeriod->size+i]=pnl_vect_get(B,i);
	}
	for(i=0;i<s->m;i++)
       for(j=0;j<s->n;j++)
	   {
	      optu_params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+numPeriod->size+B->size+s->n*i+j]=pnl_mat_get(s,i,j);
	   }
	for(i=0;i<c->m;i++)
       for(j=0;j<c->n;j++)
	   {
	      optu_params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+numPeriod->size+B->size+s->n*s->m+c->n*i+j]=pnl_mat_get(c,i,j);
	   }
	for(i=0;i<Vcoeff_dz->m;i++)
       for(j=0;j<Vcoeff_dz->n;j++)
	   {
	      optu_params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+numPeriod->size+B->size+s->n*s->m+c->n*c->m+Vcoeff_dz->n*i+j]=pnl_mat_get(Vcoeff_dz,i,j);
	      optu_params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+numPeriod->size+B->size+s->n*s->m+c->n*c->m+Vcoeff_dz->n*Vcoeff_dz->m+i*Vcoeff_dz->n+j]=pnl_mat_get(Vcoeff_dz1,i,j);
	   }
	for(i=0;i<U->size;i++)
    {
	   optu_params[n+attpts->size+ub->n*ub->m+ub->n*ub->m+numPeriod->size+B->size+s->n*s->m+c->n*c->m+Vcoeff_dz->n*Vcoeff_dz->m+Vcoeff_dz->n*Vcoeff_dz->m+i]=pnl_vect_get(U,i);
	}
    optu_func->params=(void *)optu_params;
	lower_bounds=pnl_vect_create(mu_vect->size);
	upper_bounds=pnl_vect_create(mu_vect->size);
	print_inner_steps=0;

	sum=0;
    for(i=0;i<ub->m;i++)
	   for(j=0;j<ub->n;j++)
		sum+=fabs(pnl_mat_get(ub,i,j))+fabs(pnl_mat_get(lb,i,j));

	if (sum== 0)//no bounds on MtM
    {
		if (muub < 1e6)
        {
			pnl_vect_set_double(lower_bounds,-muub);
			pnl_vect_set_double(upper_bounds,muub);
            result_optu=pnl_optim_intpoints_bfgs_solve(optu_func,optu_func_grad,nl_constraints,lower_bounds,upper_bounds,mu_vect0,tolerance,
							iter_max,print_inner_steps,mu_vect);
        }
		else
        {
			pnl_vect_set_double(lower_bounds,-1e6);
			pnl_vect_set_double(upper_bounds,1e6);
            result_optu=pnl_optim_intpoints_bfgs_solve(optu_func,optu_func_grad,nl_constraints,lower_bounds,upper_bounds,mu_vect0,tolerance,
							iter_max,print_inner_steps,mu_vect);
		}
    }
	else//bounds on MtM
    {
		if (muub < 1e6)
		{
			pnl_vect_set_double(lower_bounds,0);
			pnl_vect_set_double(upper_bounds,muub);
            result_optu=pnl_optim_intpoints_bfgs_solve(optu_func,optu_func_grad,nl_constraints,lower_bounds,upper_bounds,mu_vect0,tolerance,
							iter_max,print_inner_steps,mu_vect);
		}
		else
		{
			pnl_vect_set_double(lower_bounds,0);
			pnl_vect_set_double(upper_bounds,1e6);
            result_optu=pnl_optim_intpoints_bfgs_solve(optu_func,optu_func_grad,nl_constraints,lower_bounds,upper_bounds,mu_vect0,tolerance,
							iter_max,print_inner_steps,mu_vect);
		}
	}
	for(j=0;j<mu_vect->size;j++)
	{
		pnl_vect_set(mu_vect,j,-pnl_vect_get(mu_vect,j));
	}
//optimal default intensity
/*%%%%%%%%%%%%%%%%%%%%%%%%% compute lambda %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We will compute the default intensity surface from given mu with
% step-size dt<=dz.*/
    for(j=0;j<mu0->n;j++)
	   pnl_mat_set(mu0,0,j,pnl_vect_get(mu_vect,j));
    for(j=0;j<mu0->n;j++)
  	   pnl_mat_set(mu0,1,j,pnl_vect_get(mu_vect,mu0->n+j)/pnl_vect_get(attpts,0));
	for(i=2;i<mu0->m;i++)
           for(j=0;j<mu0->n;j++)
  		pnl_mat_set(mu0,i,j,pnl_vect_get(mu_vect,i*mu0->n+j)/(pnl_vect_get(attpts,i-1)-pnl_vect_get(attpts,i-2)));
	muu=pnl_mat_create_from_double(mu0->m,mu0->n,0.0);
	mud=pnl_mat_create_from_double(mu0->m,mu0->n,0.0);
	mu=pnl_mat_create_from_double(mu0->m,mu0->n,0.0);
	if (sum> 0)
	{
       for(i=0;i<mu0->m;i++)
	     for(j=0;j<mu0->n;j++)
	     {
			pnl_mat_set(mu0,i,j,pnl_mat_get(mu0,i,j)*(pnl_mat_get(s,i,j)>0)); //set mu0=0 for s<0
	     }
       for(i=0;i<muu->m;i++)
	     for(j=0;j<muu->n && j<numMat;j++)
	     {
			pnl_mat_set(muu,i,j,pnl_mat_get(mu0,i,j));
	     }
       for(i=0;i<mud->m;i++)
	     for(j=0;j<mud->n && j+numMat<mud->n;j++)
	     {
			pnl_mat_set(mud,i,j,pnl_mat_get(mu0,i,j+numMat));
	     }
	   pnl_mat_clone(mu,mud);
       pnl_mat_minus_mat(mu,muu);
	}
	else
	{
	  pnl_mat_clone(muu,mu0);
	  pnl_mat_clone(mud,mu0);
      for(i=0;i<mu->m;i++)
	     for(j=0;j<mu->n;j++)
	     {
	    	pnl_mat_set(mu,i,j,pnl_mat_get(mu0,i,j)*(pnl_mat_get(s,i,j)>0)); //set mu=0 for s<0
	     }
	}

	Phicoeff1 = pnl_mat_create_from_double(numTran,numMat,0.0);                                                  //coefficients of Phi at t1
	Phicoeff = pnl_mat_create_from_double(numTran,numMat,0.0);                                                   //coefficients of Phi at j
	PhicoeffJ = pnl_mat_create_from_double(numTran,numMat,0.0);                                                  //coefficients of Phi at T

	for(i=0;i<numTran-1;i++)
          for(j=0;j<numMat;j++)
             pnl_mat_set(Phicoeff1,i,j,pnl_mat_get(mu,i+1,j)*(1+pnl_mat_get(s,i+1,j)*dz1-Bdz)-pnl_mat_get(mu,i+2,j)*(1+pnl_mat_get(s,i+2,j)*dz1-Bdz));
	for(i=0;i<numTran-1;i++)
          for(j=0;j<numMat;j++)
             pnl_mat_set(Phicoeff,i,j,pnl_mat_get(mu,i+1,j)*(1+pnl_mat_get(s,i+1,j)*dz-Bdz)-pnl_mat_get(mu,i+2,j)*(1+pnl_mat_get(s,i+2,j)*dz-Bdz));
	for(i=0;i<numTran-1;i++)
          for(j=0;j<numMat;j++)
             pnl_mat_set(PhicoeffJ,i,j,pnl_mat_get(mu,i+1,j)*(1+pnl_mat_get(s,i+1,j)*dz)-pnl_mat_get(mu,i+2,j)*(1+pnl_mat_get(s,i+2,j)*dz));

        for(j=0;j<numMat;j++)
           pnl_mat_set(Phicoeff1,numTran-1,j,pnl_mat_get(mu,0,j)*(pnl_mat_get(s,0,j)*dz1/(1-recov)+1-Bdz)+pnl_mat_get(mu,mu->m-1,j)*(1+pnl_mat_get(s,s->m-1,j)*dz1-Bdz));
        for(j=0;j<numMat;j++)
           pnl_mat_set(Phicoeff,numTran-1,j,pnl_mat_get(mu,0,j)*(pnl_mat_get(s,0,j)*dz/(1-recov)+1-Bdz)+pnl_mat_get(mu,mu->m-1,j)*(1+pnl_mat_get(s,s->m-1,j)*dz-Bdz));
        for(j=0;j<numMat;j++)
           pnl_mat_set(PhicoeffJ,numTran-1,j,pnl_mat_get(mu,0,j)*(pnl_mat_get(s,0,j)*dz/(1-recov)+1)+pnl_mat_get(mu,mu->m-1,j)*(1+pnl_mat_get(s,s->m-1,j)*dz));

	Phi=pnl_mat_create_from_double(pnl_iround(pnl_vect_get(numPeriod,numPeriod->size-1)+1)*numMat,(numName+1),0.0); //compute (negative) Phi for all maturities

	temp_vect1=pnl_vect_create(PhicoeffJ->m);
	temp_vect2=pnl_vect_create(numName+1);
	temp_mat=pnl_mat_create(c->n,c->m);
	pnl_mat_tr(temp_mat,c);
	iMat = numMat-1;
	while (iMat >= 0)
	{
		pnl_mat_get_col(temp_vect1,PhicoeffJ,iMat);
		temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);

		for(j=0;j<Phi->n;j++)
			pnl_mat_set(Phi,iMat*(Phi->m/numMat)+pnl_iround(pnl_vect_get(numPeriod,iMat)),j,-pnl_vect_get(B,pnl_iround(pnl_vect_get(numPeriod,iMat))-1)
			*pnl_vect_get(temp_vect2,j));

		pnl_mat_get_col(temp_vect1,Phicoeff1,iMat);
		temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);

		for (i =pnl_iround(pnl_vect_get(numPeriod,iMat))-1;i>1;i--)
			for(j=0;j<Phi->n;j++)
				pnl_mat_set(Phi,iMat*(Phi->m/numMat)+i,j,-pnl_vect_get(B,i-1)*pnl_vect_get(temp_vect2,j));

		pnl_mat_get_col(temp_vect1,Phicoeff1,iMat);
		temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);

		for(j=0;j<Phi->n;j++)
			pnl_mat_set(Phi,iMat*(Phi->m/numMat)+1,j,-pnl_vect_get(B,0)*pnl_vect_get(temp_vect2,j));

		sum1=0;
		for(i=0;i<mu->m-2;i++)
			sum1+=pnl_mat_get(mu,i+2,iMat)*(pnl_vect_get(attpts,i+1)-pnl_vect_get(attpts,i));

		sum2=0;
		for(i=1;i<pnl_iround(pnl_vect_get(numPeriod,iMat));i++)
			sum2+=pnl_vect_get(B,i);

		for(j=0;j<Phi->n;j++)
			pnl_mat_set(Phi,iMat*(Phi->m/numMat),j,-(-pnl_vect_get(B,0)*(sum1+pnl_mat_get(mu,1,iMat)*pnl_vect_get(attpts,0)+pnl_mat_get(mu,0,iMat))
			+pnl_mat_get(mu,1,iMat)*pnl_vect_get(U,iMat)*pnl_vect_get(attpts,0)-pnl_mat_get(mu,0,iMat)*pnl_mat_get(s,0,iMat)*recov/(1-recov)
			*(pnl_vect_get(B,0)*dz1+sum2)*dz));

		iMat = iMat - 1;
	}
	sum=0;
    for(i=0;i<ub->m;i++)
	   for(j=0;j<ub->n;j++)
		sum+=pnl_mat_get(ub,i,j)*pnl_mat_get(muu,i,j)-pnl_mat_get(lb,i,j)*pnl_mat_get(mud,i,j);
	for(i=1;i<=numMat;i++)
		for(j=0;j<Phi->n;j++)
			pnl_mat_set(Phi,i*(Phi->m/numMat)-1,j,pnl_mat_get(Phi,i*(Phi->m/numMat)-1,j)+pnl_mat_get(Phi,(i-1)*(Phi->m/numMat),j)+sum);

	ePhi=pnl_mat_create(Phi->m/numMat,numName+1);
    for(i=0;i<ePhi->m;i++)
		for(j=0;j<ePhi->n;j++)
		{
			sum=0;
			for(k=0;k<numMat;k++)
				sum+=pnl_mat_get(Phi,k*(Phi->m/numMat)+i,j);
			pnl_mat_set(ePhi,i,j,exp(sum));//exp of sum of Phi in different maturities
		}
	u=pnl_mat_create(pnl_iround(pnl_vect_get(T,T->size-1)/dt)+1,ePhi->n);
	temp_vect1=pnl_vect_create(u->n);
	temp_vect2=pnl_vect_create(u->n);
	temp_mat=pnl_mat_create(Vcoeff_dt->n,Vcoeff_dt->m);
	pnl_mat_tr(temp_mat,Vcoeff_dt);

	pnl_mat_get_row(temp_vect1,ePhi,ePhi->m-1);
	pnl_mat_set_row(u,temp_vect1,u->m-1);

	temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);
	pnl_mat_set_row(u,temp_vect2,u->m-2);
	for (i=u->m-3;i>=pnl_iround(trunc(dz1/dt));i--)
	{
		if(((i+1)*dt-dz1)-dz*pnl_iround(trunc(((i+1)*dt-dz1)/dz))==0)//terminal value on payment dates
		{
			for(j=0;j<u->n;j++)
				pnl_vect_set(temp_vect1,j,(pnl_mat_get(u,i+1,j)*pnl_mat_get(ePhi,pnl_iround((((i+1)*dt-dz1)/dz)+1),j)));
			temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);
			pnl_mat_set_row(u,temp_vect2,i);
		}
		else//terminal value on intermediate dates
		{
			for(j=0;j<u->n;j++)
				pnl_vect_set(temp_vect1,j,pnl_mat_get(u,i+1,j));
			temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);
			pnl_mat_set_row(u,temp_vect2,i);
		}
	}

	for(j=0;j<u->n;j++)
		pnl_vect_set(temp_vect1,j,(pnl_mat_get(u,pnl_iround((dz1/dt)),j)*pnl_mat_get(ePhi,1,j)));
	temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);
	pnl_mat_set_row(u,temp_vect2,pnl_iround((dz1/dt))-1);
	for (i=pnl_iround((dz1/dt))-2;i>=0;i--)
	{
		for(j=0;j<u->n;j++)
			pnl_vect_set(temp_vect1,j,pnl_mat_get(u,i+1,j));
		temp_vect2=pnl_mat_mult_vect(temp_mat,temp_vect1);
		pnl_mat_set_row(u,temp_vect2,i);
	}
	//portfolio default intensity
	for(i=0;i<lambda->m;i++)
		pnl_mat_set(lambda,i,numName,0.0);
	for(i=0;i<lambda->m;i++)
	   for(j=0;j<numName;j++)
		pnl_mat_set(lambda,i,j,pnl_vect_get(lam0,j)*pnl_mat_get(u,i,j+1)/pnl_mat_get(u,i,j));

	pnl_vect_free(&numPeriod);
	pnl_vect_free(&T);
	pnl_vect_free(&B);
	pnl_vect_free(&U);
	pnl_vect_free(&lower_bounds);
	pnl_vect_free(&upper_bounds);
	pnl_mat_free(&mu);
	pnl_mat_free(&mu0);
	pnl_vect_free(&mu_vect);
	pnl_vect_free(&mu_vect0);
	pnl_mat_free(&c);
	pnl_mat_free(&Vcoeff_dz);
	pnl_mat_free(&Vcoeff_dz1);
	pnl_mat_free(&Vcoeff_dt);
	pnl_mat_free(&temp1);
	pnl_mat_free(&temp2);
	pnl_mat_free(&s);
	free(nl_constraints);
    free(optu_func);
	pnl_mat_free(&mud);
	pnl_mat_free(&muu);
	pnl_mat_free(&Phicoeff1);
	pnl_mat_free(&Phicoeff);
	pnl_mat_free(&PhicoeffJ);
	pnl_mat_free(&Phi);
	pnl_mat_free(&temp_mat);
	pnl_mat_free(&ePhi);
	pnl_mat_free(&u);
	pnl_vect_free(&temp_vect1);
	pnl_vect_free(&temp_vect2);
	free(optu_params);

	return 0;
}

int main()
{
	time_t t0, tf;

	int j,
	numName = 125;//total No. of names
	double
	r = 0.03,//interest rate
	dz = 0.25,//payment time interval
	dt = 0.0005,//dynamic prog time step
    recov = 0.4;//recovery rate
	//Optimization parameters
	double
	tolerance=1e-6,
    iter_max=1000,
	muub=1e6;

    PnlMat *mktsprd,*ub,*lb,*mu0,*lambda;
    PnlVect *T,*attpts,*lam0;

	srand(time(0));
    t0 = time (NULL);

	mktsprd=pnl_mat_create(7,1);
    pnl_mat_set(mktsprd,0,0,-999);pnl_mat_set(mktsprd,1,0,29.875);pnl_mat_set(mktsprd,2,0,98);pnl_mat_set(mktsprd,3,0,34.5);
    pnl_mat_set(mktsprd,4,0,14);pnl_mat_set(mktsprd,5,0,8.125);pnl_mat_set(mktsprd,6,0,3.125);

    T=pnl_vect_create(1);
    pnl_vect_set(T,0,5);

    attpts=pnl_vect_create(6);
    pnl_vect_set(attpts,0,0.03);pnl_vect_set(attpts,1,0.06);pnl_vect_set(attpts,2,0.09);pnl_vect_set(attpts,3,0.12);
	pnl_vect_set(attpts,4,0.22);pnl_vect_set(attpts,5,1);

    ub=pnl_mat_create_from_double(7,1,0.0);
	lb=pnl_mat_create_from_double(7,1,0.0);
	mu0=pnl_mat_create_from_double(7,1,0.0);

    lam0 = pnl_vect_create_from_double(numName+1,1.0);
	pnl_vect_set(lam0,numName,0.);

    lambda=pnl_mat_create(pnl_iround(pnl_vect_get(T,T->size-1)/dt)+1,numName+1);

    calIntensity(r,numName,dz,dt,recov,T,mktsprd,attpts,lam0,ub,lb,mu0,muub,tolerance,iter_max,lambda);

	for(j=0;j<lambda->n;j++)
	{
		printf("%f  ",pnl_mat_get(lambda,49,j));
	}
	printf("\n\n\n");

    pnl_mat_resize(mktsprd,7,3);

    pnl_mat_set(mktsprd,0,0,-999);pnl_mat_set(mktsprd,1,0,29.875);pnl_mat_set(mktsprd,2,0,98);pnl_mat_set(mktsprd,3,0,34.5);
    pnl_mat_set(mktsprd,0,1,-999);pnl_mat_set(mktsprd,1,1,47.55);pnl_mat_set(mktsprd,2,1,196.5);pnl_mat_set(mktsprd,3,1,54.5);
    pnl_mat_set(mktsprd,0,2,-999);pnl_mat_set(mktsprd,1,2,58.75);pnl_mat_set(mktsprd,2,2,512.5);pnl_mat_set(mktsprd,3,2,103);
    pnl_mat_set(mktsprd,4,0,14);pnl_mat_set(mktsprd,5,0,8.125);pnl_mat_set(mktsprd,6,0,3.125);
    pnl_mat_set(mktsprd,4,1,31.5);pnl_mat_set(mktsprd,5,1,13.5);pnl_mat_set(mktsprd,6,1,6.25);
    pnl_mat_set(mktsprd,4,2,51.5);pnl_mat_set(mktsprd,5,2, 23.55);pnl_mat_set(mktsprd,6,2,9.5);

    pnl_vect_resize(T,3);
    pnl_vect_set(T,0,5);pnl_vect_set(T,1,7);pnl_vect_set(T,2,10);

    pnl_vect_resize(attpts,6);
    pnl_vect_set(attpts,0,0.03);pnl_vect_set(attpts,1,0.06);pnl_vect_set(attpts,2,0.09);pnl_vect_set(attpts,3,0.12);pnl_vect_set(attpts,4,0.22);pnl_vect_set(attpts,5,1);

	pnl_mat_resize(ub,7,3);
	pnl_mat_resize(lb,7,3);
    pnl_mat_set_double(ub,0.0);
	pnl_mat_set_double(lb,0.0);
	pnl_mat_resize(mu0,7,3);
	pnl_mat_set_double(mu0,0.0);

    pnl_mat_get_row(lam0,lambda,3999);
    pnl_mat_resize(lambda,pnl_iround(pnl_vect_get(T,T->size-1)/dt)+1,numName+1);

    calIntensity(r,numName,dz,dt,recov,T,mktsprd,attpts,lam0,ub,lb,mu0,muub,tolerance,iter_max,lambda);
	for(j=0;j<lambda->n;j++)
	{
		printf("%f  ",pnl_mat_get(lambda,49,j));
	}
	printf("\n\n");
	tf=time (NULL);
	printf("%d\n",(int)(tf-t0));

	pnl_vect_free(&attpts);
	pnl_mat_free(&ub);
	pnl_mat_free(&lb);
	pnl_mat_free(&lambda);
	pnl_mat_free(&mu0);
	pnl_mat_free(&mktsprd);
	pnl_vect_free(&T);
	pnl_vect_free(&lam0);

	return 0;

}


