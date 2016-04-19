#include "utils1.h"
// ===================================================================
// Loi uniforme
// ===================================================================
double uniforme(void) { return (double) rand()/RAND_MAX;}
// ===================================================================
//==================================================================================
// Affiche_Sol
//==================================================================================
void affiche_sol(int dimx,double *sol,double *x,double fx,double *grad,double *grad0)
{
  int i;
  //
  // sol
  printf("sol  = ");
  for(i=0;i<dimx;i++) printf("%f ",sol[i]);
  printf("\n");
  printf("x*  = ");
  for(i=0;i<dimx;i++) printf("%f, ",x[i]); printf("\n");
  printf("fx = %f\n",fx);
  printf("g  = ");
  for(i=0;i<dimx;i++) printf("%f ",grad[i]);printf("\n");
  printf("g0  = ");
  for(i=0;i<dimx;i++) printf("%f ",grad0[i]);printf("\n");
  printf("g/g0  = ");
  for(i=0;i<dimx;i++) printf("%f ",fabs(grad[i]/grad0[i]));
  printf("\n");
  //
  //
}
//==================================================================================
// F_Affiche_Sol
//==================================================================================
void faffiche_sol(FILE *ftest,int dimx,double *sol,double *x,double fx,double *grad,double *grad0)
{
  int i;
  //
  // sol
  fprintf(ftest,"sol  = ");
  for(i=0;i<dimx;i++) fprintf(ftest,"%f ",sol[i]);
  fprintf(ftest,"\n");
  fprintf(ftest,"x*  = ");
  for(i=0;i<dimx;i++) fprintf(ftest,"%f, ",x[i]); fprintf(ftest,"\n");
  fprintf(ftest,"fx = %f\n",fx);
  fprintf(ftest,"g  = ");
  for(i=0;i<dimx;i++) fprintf(ftest,"%f ",grad[i]);fprintf(ftest,"\n");
  fprintf(ftest,"g/g0  = ");
  for(i=0;i<dimx;i++) fprintf(ftest,"%f ",fabs(grad[i]/grad0[i]));
  fprintf(ftest,"\n");
  //
  //
}
//==================================================================================
// Norme
//==================================================================================
double norme(int dimx,double *x)
{
  int i;
  double norm;
  //
  norm = 0.;
  //
  for(i=0;i<dimx;i++) norm+=fabs(x[i])*fabs(x[i]);
  //
  return sqrt(norm);
}
//
double dist(int dimx,double *x,double *y)
{
  int i;
  double norm;
  //
  norm = 0.;
  //
  for(i=0;i<dimx;i++) norm+=fabs(x[i]-y[i])*fabs(x[i]-y[i]);
  //
  return sqrt(norm);
}
//==================================================================================
// Init_Sol
//==================================================================================
void init_sol(int TypeModel,int *dimx,double *x,int dimx1,int dimx2,int dimx3,double V0,double kappa,double theta,double sigmav,double rho,double lambda, double m0,double v)
{
  
  if(TypeModel==1 || TypeModel==111)
	{
	  *dimx = dimx1;
	  x[0] = V0;
	  x[1] = kappa;
	  x[2] = theta;
	  x[3] = sigmav;
	  x[4] = rho;
	}
  else if(TypeModel==11)
	{
	  *dimx = dimx1;
	  x[0] = V0;
	  x[1] = kappa/sigmav;
	  x[2] = theta;
	  x[3] = sigmav;
	  x[4] = rho;
	}
  else if(TypeModel==2)
	{
	  *dimx = dimx2;
	  x[0] = V0;
	  x[1] = lambda;
	  x[2] = m0;
	  x[3] = v;
	}
  else if(TypeModel==3)
	{
	  *dimx = dimx3;
	  x[0] = V0;
	  x[1] = kappa;
	  x[2] = theta;
	  x[3] = sigmav;
	  x[4] = rho;
	  x[5] = lambda;
	  x[6] = m0;
	  x[7] = v;
	}
  else
	{
	  *dimx = dimx3;
	  x[0] = V0;
	  x[1] = kappa/sigmav;
	  x[2] = theta;
	  x[3] = sigmav;
	  x[4] = rho;
	  x[5] = lambda;
	  x[6] = m0;
	  x[7] = v;
	}
  //  
}
//==================================================================================
// Init_Bornes
//==================================================================================
void init_bornes(int TypeModel,int dimx,int *nbd,double *xmin,double *xmax)
{
  int i;
  
  for(i=0;i<dimx;i++) {nbd[i] = 2;}
  if(TypeModel==1 || TypeModel==11 || TypeModel==111)
	{
	  xmin[0] = 25.e-4;  xmax[0] = 1.;  nbd[0] =2;
	  //	   	  	  xmin[0] = V0;  xmax[0] = V0;
	  xmin[1] = 0.01;  xmax[1] = 5.;  nbd[1] = 2;
	  xmin[2] = 1.e-3;  xmax[2] =1.; nbd[2] = 2;
	  xmin[3] = 0.05;  xmax[3] = 1.; nbd[3] = 2;
	  xmin[4] = -1.;    xmax[4] = 1.;
	}
  else if(TypeModel==2)
	{
	  xmin[0] = 25.e-4;  xmax[0] = 1.;
	  //  xmin[0] = 0.0005;  xmax[0] = 0.06;
	  xmin[1] = 5.e-2;   xmax[1] = 2.; nbd[1] = 2;
	  //  xmin[1] =5.e-1;   xmax[1] = 1.5; nbd[1] = 1;
	  //	  xmin[2] = -.9;    xmax[2] = .4; nbd[2] = 0;
	  xmin[2] = -1.;    xmax[2] = 1.; nbd[2] = 1;
	  xmin[3] = 1.e-2;  xmax[3] = 2.; nbd[3] = 1;
	}
  else
	{
	  xmin[0] = 25.e-4;  xmax[0] = .5;
	  xmin[1] = 1.e-3;  xmax[1] = 10.;   nbd[1] = 2;
	  xmin[2] = 1.e-3;  xmax[2] = 1.;  nbd[2] = 2;
	  xmin[3] = 1.e-3;  xmax[3] = 1.;
	  xmin[4] = -1.;    xmax[4] = 1.;
	  xmin[5] = 1.e-2;   xmax[5] = 10.; nbd[5] = 1;
	  xmin[6] = -2.;    xmax[6] = 2.;   nbd[6] = 1;
	  xmin[7] = 1.e-2;  xmax[7] = 2.;   nbd[7] = 1; 
	}
}
//==================================================================================
// Init_Prix_Obs
//==================================================================================
void init_PrixObs(int nbdata1,int nbdata2,int TypeModel,double Kmin,double Kmax,double Tmin,double Tmax,double St0,double r,double divid,double V0,double kappa,double theta,double sigmav,double rho,double lambda, double m0,double v)
{
  int ibidon,i,ii,j,j2,type; 
  double K,dK,T,dT,wi=1.;
  double prix,delta,sigma_imp_obs;
  NumFunc_1  *p;
  //
  FILE *fich;
  //
  fich = fopen ("DataMarketSynthetiques.txt", "wt");
  // 
  p=(NumFunc_1*)malloc(sizeof(NumFunc_1));
  //
  dK   = (Kmax-Kmin)/(nbdata2-1);
  dT   = (Tmax-Tmin)/(nbdata1-1);
  //
  for(j=0;j<nbdata1;j++)
    {
      //
      T = Tmin + j*dT;
      K = Kmin - dK;
      //
      for(i=0;i<nbdata2;i++)
	{
	  ii = j*nbdata2+i;
	  // on fait varier K;
	  K               = K + dK;
	  //		  type = j%2+1;   // 1 pour call et 2 pour put
	  //
	  /*
	    if(K <= St0)
	    type =  1 ;
	  else
	    type =  2 ;
	  */
	  type =1;
	  //
	  p->Par[0].Val.V_DOUBLE =K;
	  //
	  if(type==1) 
	    {
	      //
	      if(TypeModel==1)
		{
		  ibidon=FT_Call_Heston(St0,p,T,r,divid,V0,kappa,theta,sigmav,rho,&prix,&delta);
		}
	      else if(type==2)
		{
		  ibidon=FT_Call_Merton(St0,p,T,r,divid,V0,lambda,m0,v,&prix,&delta);
		}
	      else
		{
		  ibidon=FT_Call_HestMert(St0,p,T,r,divid,V0,kappa,theta,sigmav,rho,lambda,m0,v,&prix,&delta);
		}
	      //
	    }
	  else 
	    {
	      //
	      if(TypeModel==1)
		{
		  ibidon=FT_Put_Heston(St0,p,T,r,divid,V0,kappa,theta,sigmav,rho,&prix,&delta);
		}
	      else if(TypeModel==2)
		{
		  ibidon=FT_Put_Merton(St0,p,T,r,divid,V0,lambda,m0,v,&prix,&delta);
		}
	      else
		{
		  ibidon=FT_Put_HestMert(St0,p,T,r,divid,V0,kappa,theta,sigmav,rho,lambda,m0,v,&prix,&delta);
		}
	      //
	    }
	  //		  
	  // calcul de la vol_implicite
	  sigma_imp_obs = SigmaImplicite(eps,a0,b0,type,St0,T,K,r,divid,prix,&j2);
	  if(j2 == -1) {
	    printf("Warring : Les prix calcule par le modele ne sont pas compatibles avec BS. \n"); 
	    printf("          Pour St0=%f,T=%f,K=%f,sigma_implicite=%f \n",St0,T,K,sigma_imp_obs);
	  }
	  //
	  fprintf(fich,"%d %f %f %f %f %f %f \n",type,T,K/100.,sigma_imp_obs,r,divid,wi);
	  //
	}
    }
  //
  free(p);
  //
  fclose(fich);
  //
}
//==================================================================================
// Trace_1d
//==================================================================================
void trace_1d(int iplot1,int iplot2,int nbplot,int dimx,double *x0,double *y0,double *y1,double *xmin,double *xmax)
{
  int iplot,i0,i,ii;
  double x[8];
  FILE *fichier;
  double dx,dy,x10,x20;
  double pt1x,pt1y,pt2x,pt2y,pta;

  //
  for(i=0;i<dimx;i++) x[i] = x0[i];
  //
  pt1x = 0.024600;
  pt1y = 0.5;
  pt2x = 0.001000;
  pt2y = 1.5;
  pta  = (pt2y-pt1y)/(pt2x-pt1x);
  //
  iplot  = iplot1; 
  i0 = 0;
  
  // dx = i0*(xmax[iplot] - xmin[iplot])/ (2.*nbplot) + (1-i0)*(xmax[iplot] - xmin[iplot])/nbplot;
  //x0 = x[iplot];
  
  dx = i0*(xmax[iplot1] - xmin[iplot1])/ (2.*nbplot) + (1-i0)*(xmax[iplot1] - xmin[iplot1])/nbplot;
  dy = i0*(xmax[iplot2] - xmin[iplot2])/ (2.*nbplot) + (1-i0)*(xmax[iplot2] - xmin[iplot2])/nbplot;
  x10 = x[iplot1];
  x20 = x[iplot2];
  printf("x10=%f, x20=%f\n",x10,x20) ;

  //
  //  for(i=0;i<dimx;i++) x[i] = 2.*x[i];
  
  //
  fichier = fopen ("prix.dat", "wt");
  for(i=-nbplot*i0;i<=nbplot;i++)
		{
 		  x[iplot]    = (x10 + (double)i*dx)*i0 +  (xmin[iplot] + (double)i*dx)*(1-i0);
 		  x[iplot2]    = pt1y + pta*(x[iplot]-pt1x);
		  ii          = nbplot*i0+i;
		  y0[ii]      = x[iplot];
		  y1[ii]      = costFunction(dimx,x);
		  //		  fprintf(fichier,"%f %f \n",y0[ii],y1[ii]);
		  fprintf(fichier,"%f %f %f \n",y0[ii],x[iplot2],y1[ii]);
		  if(x[iplot2]<0.9*xmin[iplot2] || x[iplot2]>1.1*xmax[iplot2]) break;
		  
	  }
  fclose(fichier);

  
 
  //plot(nbplot,y0,y1);
}
//==================================================================================
// Trace_3d
//==================================================================================
void trace_3d(int iplot1,int iplot2,int nbplot,int nbdata1,int dimx,double *x0,double *y0,double *y1,double fy0y1[2*nbplot+1][2*nbplot+1],double *xmin,double *xmax)
{
  int i0,i,ii,j,jj;
  double *fxy;
  double dx,dy,x10,x20;
  double x[8];
  FILE *fichier;
  //
  for(i=0;i<dimx;i++) x[i] = x0[i];
  //
  fxy = (double *) malloc((2*nbdata1+1)*(2*nbdata1+1)*sizeof(double));

  i0 = 0;
  
   
  dx = i0*(xmax[iplot1] - xmin[iplot1])/ (2.*nbplot) + (1-i0)*(xmax[iplot1] - xmin[iplot1])/nbplot;
  dy = i0*(xmax[iplot2] - xmin[iplot2])/ (2.*nbplot) + (1-i0)*(xmax[iplot2] - xmin[iplot2])/nbplot;
  x10 = x[iplot1];
  x20 = x[iplot2];

	 
  for(i=-nbplot*i0;i<=nbplot;i++)
	for(j=-nbplot*i0;j<=nbplot;j++)
	  {
		{
		  //x[iplot1]   = xmin[iplot1] + (double)i*dx;
		  //x[iplot2]   = xmin[iplot2] + (double)j*dy;
		  x[iplot1]   = (x10 + (double)i*dx)*i0 +  (xmin[iplot1] + (double)i*dx)*(1-i0);
		  x[iplot2]   = (x20 + (double)j*dy)*i0 +  (xmin[iplot2] + (double)j*dy)*(1-i0);
		  //
		  ii = nbplot*i0+i;
		  jj = nbplot*i0+j;
		  //
		  y0[ii]       = x[iplot1];
		  y1[jj]       = x[iplot2];
		  fy0y1[ii][jj] = costFunction(dimx,x);
		  fxy[jj*(nbplot*i0+nbplot+1)+ii+1] = fy0y1[ii][jj];
		  //printf("x10=%f, x20=%f, fx=%f \n",x[iplot1],x[iplot2],fy0y1[ii][jj]) ;	  
		}
	  }
  
  //  plot(nbplot,y0,y1);
  //
  fichier = fopen ("surfprix.dat", "wt");
  //  fich = fopen ("logsurfprix.dat", "wt");
  //  for(i=0;i<nbplot;i++)
  //	for(j=0;j<nbplot;j++)
  for(i=-nbplot*i0;i<=nbplot;i++)
	for(j=-nbplot*i0;j<=nbplot;j++)
	  {
		{
		  //
		  ii = nbplot*i0+i;
		  jj = nbplot*i0+j;
		  //		  
		  fprintf(fichier,"%f %f %f \n",y0[ii],y1[jj],fy0y1[ii][jj]);
		}
	  }
  fclose(fichier);
  //fclose(fich);
  
  //  plot3d(nbplot,nbplot,y0,y1,fy0y1);
  //  myplot3d(nbplot+1,nbplot+1,y0,y1,fxy);

  free(fxy);
  
}
//==================================================================================
// Convert_date
//==================================================================================
void convert_date(int jj,int mm, int aa,char *date)
{
  //
  aa = aa + 1900;
  if(jj<10 && mm < 10 )
    sprintf(date,"0%d_0%d_%d",jj,mm,aa);
  else if(jj<10 && mm >= 10 )
	sprintf(date,"0%d_%d_%d",jj,mm,aa);
  else if(jj>=10 && mm < 10 )
	sprintf(date,"%d_0%d_%d",jj,mm,aa);
  else 
    sprintf(date,"%d_%d_%d",jj,mm,aa);
  //
}
//==================================================================================
//
//==================================================================================

void readDataMarket(DataMarket *DM)
{
  int i,j,type;
  double T,K,St0,r,d,wi,sigma,prix,delta;
  FILE *fichier;
  //
  fichier = fopen ("DataMarket.txt", "r");
  //
  for(i=0;i<DM->nbdata;i++)
	{
	  //
	  fscanf(fichier,"%d %lf %lf %lf %lf %lf %lf\n",&type,&T,&K,&sigma,&r,&d,&wi);
	  //
	  r = log(1.+r);
	  d = log(1.+d/100.);
	  //
	  St0 = 100.;
	  //
	  DM->TypeOpt[i]   = type;
	  DM->St0[i] = St0;
	  DM->T[i]   = T;
	  DM->K[i]   = St0*K;
	  DM->r[i]   = r;
	  DM->d[i]   = d;
	  DM->wi[i]   = wi;
	  DM->SigmaImpObs[i]   = sigma;
	  //
	  j =  MyCall_BlackScholes_73(St0,St0*K,T,r,d,sigma,&prix,&delta);
	  DM->PrixObs[i]   = prix;
	  //
	}
  //
  fclose(fichier);
  //
}
//==================================================================================                                                      
//                                                                                                                                        
//==================================================================================
void read_params(int TypeModel,int *typen,int *nbetapes,int *nbdata,double *xinit)
{
  //
  FILE *fichier;
  //
  fichier = fopen ("in.dat","r");
  //
  fscanf(fichier,"%d \n",typen);
  fscanf(fichier,"%d \n",nbetapes);
  fscanf(fichier,"%d \n",nbdata);
  if(TypeModel==1)
    {
      fscanf(fichier,"%lf %lf %lf %lf %lf\n",&xinit[0],&xinit[1],&xinit[2],&xinit[3],&xinit[4]);
    }
  else if(TypeModel==2)
    {
      fscanf(fichier,"%lf %lf %lf %lf\n",&xinit[0],&xinit[1],&xinit[2],&xinit[3]);
    }
  else
    {
      fscanf(fichier,"%lf %lf %lf %lf %lf %lf %lf %lf\n",&xinit[0],&xinit[1],&xinit[2],&xinit[3],&xinit[4],&xinit[5],&xinit[6],&xinit[7]);
    }
  //  printf("%f %f %f %f %f\n",xinit[0],xinit[1],xinit[2],xinit[3],xinit[4]);
  //
  fclose(fichier);
  //
}
//==================================================================================                                                      
//                                                                                                                                        
//==================================================================================
void ReInitParams(int nbdata,DataMarket *DM,int TypeNorme,int TypeModel,int LogNorme,int MelangeNormes)
{
  //
  FreeCostFunction();
  FreeGradFunction();
  initParams(nbdata,DM,TypeNorme,TypeModel,LogNorme,MelangeNormes);
  initParamsGrad(nbdata,DM,TypeNorme,TypeModel,LogNorme,MelangeNormes);
  //
}
//==================================================================================                                                      
//                                                                                                                                        
//==================================================================================
void BloqueVariables(int dimx,int TypeModel,int lequel,int *nbd,double *x,double *xmin,double *xmax)
{
  //
  int i;
  //
  switch(TypeModel)
    {
      // Heston
    case 1:
      //
      printf(" On bloque les variables V0, theta et rho \n");
      for(i=0;i<dimx;i++) {
	if(i!=1 && i!=3) {
	  xmin[i] = x[i];
	  xmax[i] = xmin[i];
	  nbd[i]  = 2;
	}
      }
      //    
      break;
      // Merton
    case 2:
      //
      printf(" On bloque les variables V0 et Lambda \n");
      for(i=0;i<dimx;i++) {
	//	if(i!=2 && i!=3) {
	if(i==0 || i==1) {
	  xmin[i] = x[i];
	  xmax[i] = xmin[i];
	  nbd[i]  = 2;
	}
      }
      //
      break;
      // Merton+Heston
    case 3:
      if(lequel==1)
	{
	  printf(" On bloque les variables Heston\n");
	  for(i=0;i<dimx;i++) {
	    if(i==1 || i==2 || i==3 || i==4) {
	      xmin[i] = x[i];
	      xmax[i] = xmin[i];
	      nbd[i]  = 2;
	    }
	  }
	}
      else
	{
	  printf(" On bloque les variables Merton\n");
	  for(i=0;i<dimx;i++) {
	    if(i==5 || i==6 || i==7) {
	      xmin[i] = x[i];
	      xmax[i] = xmin[i];
	      nbd[i]  = 2;
	    }
	  }

	}
      break;
    }
  //
}
//==================================================================================                                                      
//                                                                                                                                        
//==================================================================================
void BloqueVariable_i(int dimx,int i,int *nbd,double *x,double *xmin,double *xmax)
{
  //
  int j;
  //
  printf(" On bloque toutes les variables sauf %d \n",i);
  //
  for(j=0;j<dimx;j++) 
    if(j!=i)
      {
	xmin[j] = x[j];
	xmax[j] = xmin[j];
	nbd[j]  = 2;
      }
  //
  //
}
//==================================================================================                                                      
//                                                                                                                                        
//==================================================================================
void ReInitBornes(int dimx,int *nbd0,double *xmin0,double *xmax0,int *nbd,double *xmin,double *xmax)
{
  //
  int i;
  //
  for(i=0;i<dimx;i++) {
    xmin[i] = xmin0[i];
    xmax[i] = xmax0[i];
    nbd[i]  = nbd0[i];
  }
  //
}
//==================================================================================                                                      
//                                                                                                                                        
//==================================================================================
