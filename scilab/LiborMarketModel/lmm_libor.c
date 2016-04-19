
#include"lmm_header.h"


int mallocLibor(Libor **ptLib ,int numOfMat, double tenorVal )
{
  int i;
  Libor *pt;
  pt=(Libor*)malloc(sizeof(Libor));
  
  pt->numberOfMaturities=numOfMat;
  pt->tenor=tenorVal;

  pt->libor=(double*)malloc(sizeof(double)*pt->numberOfMaturities);
  pt->maturity=(double*)malloc(sizeof(double)*pt->numberOfMaturities);
  for(i=0;i<pt->numberOfMaturities;i++)
    {
      pt->maturity[i]=i*pt->tenor;
      pt->libor[i]=0.05;
    }

  *ptLib=pt;
  return(EXIT_SUCCESS);
}


int freeLibor(Libor **ptLib)
{
  Libor  *pt;

  pt=(*ptLib);
  *ptLib=NULL;
  free(pt->libor);
  free(pt->maturity);
  
  return(1);
}

int initLibor(Libor *ptLib)
{
  int i;

  for (i=0; i<ptLib->numberOfMaturities; i++)
    {
      ptLib->libor[i]=5.0/100.;
    }
};

int readLiborFromFile(Libor **ptLib, char *fileName)
{
  /*lie les donnees des libors initiaux et leurs maturités, les donnees lues prevalent sur celles donnees dans 'initLiborList()'*/
  int i,n, etat;
  char ligne[20];
  char* pligne;
  double t, l, Tprev, delta, deltaprev;
  FILE *datas;
  double *L;
  double *T;
  
  datas=fopen(fileName, "r");
   
  if(datas==NULL)
    {
      printf("Le FICHIER N'A PU ETRE OUVERT. VERIFIER LE CHEMIN\n");
      exit(1);
    }

  n=0;
  Tprev=0;
  deltaprev=0;
  
  pligne=ligne;
  T=(double *)malloc(100*sizeof(double));
  L=(double *)malloc(100*sizeof(double));
  /* printf("OUVERTURE\n");*/
  
 
  while(1)
    {
      pligne=fgets(ligne, sizeof(ligne), datas);
      if(pligne==NULL) 
	break;
      else
	{
	  sscanf(ligne, "%lf t=%lf", &t, &l);
	  
	  T[n]=t;
	  L[n]=t;
	  
	  delta=t-Tprev;
	  Tprev=t;
	  if(delta!=deltaprev && n>0){printf("WARNING, NO CONSTANT TENOR IN LIBOR LIST!\n");}
	  deltaprev=delta;
	  
	  n++;
	} 
    }
  etat=fclose(datas);

  
  (*ptLib)->maturity=(double *)malloc(n*sizeof(double));
  (*ptLib)->libor =(double *)malloc(n*sizeof(double));

  for(i=0;i<n;i++)
    {
      (*ptLib)->maturity[i]=T[i];
      (*ptLib)->libor[i]=L[i];
    }
  (*ptLib)->numberOfMaturities=n;
  (*ptLib)->tenor=T[2]-T[1];         
  
  free(T);
  free(L);
   
  return(1); 
}



int putLiborToZero(Libor *ptLib, int index)
{
  ptLib->libor[index]=0.0;
};


int copyLibor(Libor *ptLibSrc , Libor *ptLibDest )
{
  int i ;
  
  ptLibDest->numberOfMaturities=ptLibSrc->numberOfMaturities;
  ptLibDest->tenor=ptLibSrc->tenor;

  for (i=0; i<ptLibSrc->numberOfMaturities; i++)
    {
      ptLibDest->libor[i]=ptLibSrc->libor[i];
      ptLibDest->maturity[i]=ptLibSrc->maturity[i];
    }
};


int printLibor(Libor *ptLib)
{
  int i,j;
  for(i=0;i<ptLib->numberOfMaturities;i++)
    {
      printf("%lf %lf \n",ptLib->maturity[i], ptLib->libor[i]);
      
    }
  printf("\n");
}


double computeSwapRate(Libor* ptLib, int o, int s,int m )
{
  // compute   (B(T_o,T_s)-B(T_o,T_m))/( sum_{i=s+1}^{m} \tau B(T_o,T_i) )=S(T_o,T_s,T_m) the forward swap rate 
  int j,k,l;
  double val=1.;
  double sum=0.0;
  double vald;

  for(k=s+1;k<=m;k++)
    {
      val=1.;
      for(l=o ; l<k ; l++)  
	{
	  val*=1./(1.+ ptLib->tenor*ptLib->libor[l]);
	}
      sum+=ptLib->tenor* val;
    }
  
  val=1.;
  for(k=o;k<m;k++)
    {
      val*=1./(1.+ ptLib->tenor*ptLib->libor[k]);
    }
  
  if (o!=s)
    {
      vald=1.;
      for(k=o;k<=(s-1);k++)
	{
	  vald*=1./(1.+ ptLib->tenor*ptLib->libor[k]);
	}
      return((vald-val)/sum);
    }
  else
    {
      return((1.-val)/sum);
    }
}

double computeSwapPrice(Libor* ptLib, Swaption* ptSwp,int o, int s,int m )
{
  // compute   B(T_o,T_s)-B(T_o,T_m)-K* sum_{i=s+1}^{m} \tau B(T_o,T_i)  
  //price at time T_o of a swap on T_s,....T_m
  int j,k,l;
  double val=1.;
  double sum=0.0;
  double vald;
  double price;
  
  //sum_{k=s+1}^{m} \tau B(T_o,T_k)
  for(k=s+1;k<=m;k++)
    {
      val=1.;
      //B(T_o,T_k) 
      for(l=o ; l<k ; l++)  
	{
	  val*=1./(1.+ ptLib->tenor*ptLib->libor[l]);
	}
      sum+=ptLib->tenor* val;
    }
  
  val=1.;
  //B(T_o,T_m)
  for(k=o;k<m;k++)
    {
      val*=1./(1.+ ptLib->tenor*ptLib->libor[k]);
    }
  
  if (o!=s)
    {
      vald=1.;
      //B(T_o,T_s)
      for(k=o;k<s;k++)
	{
	  vald*=1./(1.+ ptLib->tenor*ptLib->libor[k]);
	}
      price=vald-val-ptSwp->strike*sum;
    }
  else //B(T_o,T_s)=1
    {
      price=1.-val-ptSwp->strike*sum;
    }
  return(price);
}

double computeZeroCouponSum(Libor* ptLib, int o,int s,int m )
{
  // compute   sum_{i=s}^{m} \tau B(T_o,T_i) 
  int j,k,l;
  double val=1.;
  double sum=0.0;

  for(j=s;j<=m;j++)
    {
      val=1.;
      for(k=o;k<j;k++)
	{
	  val*=1./(1.+ ptLib->tenor*ptLib->libor[k]);
	}
      sum+=ptLib->tenor* val;
    }
  
  return( sum );
}

 
double computeZeroCoupon(Libor* ptLib, int o,int s)
{
  // compute B(T_o,T_s)
  if (o==s)
    {
      return(1);
    }
  else
    {
      return(computeZeroCouponSum( ptLib, o, s, s )/ptLib->tenor);
    }
}

/* ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
/* int mallocHistLibor(HistLibor **ptLib ,int numOfMat, double tenorVal, int numberOfTraj ) */
/* { */
/*   int i,l,j; */
/*   HistLibor *pt; */
/*   pt=(HistLibor*)malloc(sizeof(HistLibor)); */
  
/*   pt->numberOfMaturities=numOfMat; */
/*   pt->tenor=tenorVal; */
/*   pt->numberOfTrajectories=numberOfTraj; */

/*   pt->maturity=(double*)malloc(sizeof(double)*pt->numberOfMaturities); */
/*   pt->libor=(double***)malloc(sizeof(double **)*pt->numberOfTrajectories); */

/*   for(i=0 ; i<pt->numberOfTrajectories ; i++) */
/*     { */
/*       pt->libor[i]=(double**)malloc(sizeof(double*)*pt->numberOfMaturities); */
/*       for(j=0 ; j<pt->numberOfMaturities ; j++) */
/* 	{ */
/* 	  pt->libor[i][j]=(double*)malloc(sizeof(double)*pt->numberOfMaturities); */
/* 	  for(l=0 ; l<pt->numberOfMaturities ; l++) */
/* 	    { */
/* 	      pt->libor[i][j][l]=0.05; */
/* 	    } */
/* 	} */
/*     } */
  
/*   for(i=0 ; i<pt->numberOfMaturities ; i++) */
/*     { */
/*       pt->maturity[i]=i*pt->tenor; */
/*     } */

/*   *ptLib=pt; */
/*   return(EXIT_SUCCESS); */
/* } */

