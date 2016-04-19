/*------------------------------------*/
/*   BERMUDAN SWAPTION PRICER         */
/* see LMM_bermudan_swptns.dvi  for   */
/* the DOC.                           */
/*------------------------------------*/
/*  Nicola Moreni, Premia 2005        */
/*------------------------------------*/


#include "lmm_header.h"
#include "lmm_libor.h"
#include "lmm_numerical.h"
#include "lmm_basis.h"
#include "lmm_mathtools.h"
#include"lmm_bermudaprice.h"
#include <string.h>

static double (*Basis)(double *x, int i);
static int (*computeEvolution)(RandomGenerator* ptRand,Libor* ptLibOld,Libor* ptLibNew,Volatility* ptVol,double dt,double t);
static int check_parameters(float tenor , int numberTimeStep, int *numFac, double swaptionMat , double swapMat , double payoff_as_Regressor , long numberMCPaths , int *Regr_Basis_Dimension , char* basis_name , char* measure_name , char Explanatory , double strike);
// added function for interface

 
double  lmm_swaption_payer_bermudan_LS_pricer(float tenor , int numberTimeStep, int numFac, double swaptionMat , double swapMat ,    double payoff_as_Regressor , long numberMCPaths , int Regr_Basis_Dimension , char* basis_name , char* measure_name , char Explanatory , double strike)
{
  
  // int numberTimeStep=10;
  // int numFac=1;
  // double  swaptionMat in (years)
  // double  swapMat in  (years)
  // payoff_as_Regressor in (years) Maturity after which payoff is included in regression
  // long numberMCPaths number of monte carlo paths
  // int Regr_Basis_Dimension=4; //finite-dimensional approx. of L² 
  // char Explanatory = 'N' or 'B' or 'S' ---> Explanatory variable for regression B=Brownian, S=Nominal Swap Paying Value, N=Numeraire;
  //  int numMat ;       //WARNING: including T_0=0, maturities number is 'numMat+1'

  Volatility *ptVol;
  Libor *ptLib;
  Swaption *ptSwpt;


  double  priceVal=0.20;
  double K=strike;//strike
  int numMat;
  char Basis_Choice[10];  //See in the followings
  char Measure_Choice[10];// "    "    "    "
  double p;
  char auxstr[5];
  int I;
  double DT=(tenor/(double)numberTimeStep);
  
 
  check_parameters( tenor , numberTimeStep,  &numFac,  swaptionMat ,  swapMat ,  payoff_as_Regressor , numberMCPaths , &Regr_Basis_Dimension , basis_name , measure_name ,  Explanatory ,  strike);

  strcpy(Basis_Choice,basis_name); //"CanD" for canonical basis or "HerD" for hermite 
  strcpy(Measure_Choice,measure_name);//"Spot" or "Fwd" cfr numerical.c for evolution
 
  //Memory Allocation and Libor/Swap Structures Initialization
  numMat=(int)(swapMat/tenor);

  mallocLibor(&ptLib,numMat,tenor);
  mallocVolatility(&ptVol,numFac);
  mallocSwaption(&ptSwpt,swaptionMat,swapMat,priceVal,K,tenor);
  
  //Swaption Computation
  I=computeBermudeanSwaption(numberMCPaths,numberTimeStep,ptLib,ptSwpt,ptVol,Regr_Basis_Dimension, payoff_as_Regressor,Basis_Choice,Measure_Choice,Explanatory);
  
  p=ptSwpt->price;
  //Printing Results to STDOutput
  //printf("s=%d, e=%d, exdates=%d, payoffasregr=%d\n",s,ptLib->numberOfMaturities,numberOfExerciseDates,PayOff_As_Regressor);
  
  //initLibor(ptLib);
  //I=computeSwaptionPriceSpot(numberMCPaths,numberTimeStep,DT,ptLib,ptSwpt,ptVol);
  
  //freeing Memory
  freeSwaption(&ptSwpt);
  freeLibor(&ptLib);
  freeVolatility(&ptVol);
  
  return(p);
}


static int check_parameters(float tenor , int numberTimeStep, int *numFac, double swaptionMat , double swapMat , double payoff_as_Regressor , long numberMCPaths , int *Regr_Basis_Dimension , char* basis_name , char* measure_name , char Explanatory , double strike)
{
  double s;
  double M;

  if(swaptionMat>=swapMat)
    {
      printf(" swaption maturity must be lower than swap maturity !\n");
      exit(-1);
    }
  s=(int)(swaptionMat/tenor);
  if (fabs(swaptionMat-s*tenor)> 0.00001)
    {
      printf(" swaption maturity must be a multiple of period\n");
      exit(-1);
    }
  M=(int)(swapMat/tenor);
  if (fabs(swapMat-M*tenor)> 0.0)
    {
      printf(" swap maturity must be a multiple of period\n");
      exit(-1);
    }
  if((Explanatory!='N') && (Explanatory!='B')  && (Explanatory!='S') )
    {
       printf("Explanatory variable must be either B, N or S \n"); 
       exit(-1);
    }
  if( *numFac>2 )
    {
      printf("Number of factors too high switching to 2 \n");
      *numFac=2;
    }
  if( (*numFac==1) && (*Regr_Basis_Dimension>6))
    {
      printf("Regression basis dimension is too high swicthing to 6 \n");
      *Regr_Basis_Dimension=6;
    }
  if( (*numFac==2) && (*Regr_Basis_Dimension>21))
    {
      printf("Regression basis dimension is too high swicthing to 6 \n");
      *Regr_Basis_Dimension=6;
    }
  if( (payoff_as_Regressor<swaptionMat) || (payoff_as_Regressor>(swapMat-tenor)) )
    {
      printf("payoff_as_Regressor value is not valid, must be within [ swaptionMat  swapMat-tenor ] \n");
      printf("chosing payoff_as_Regressor=swapMat-tenor will not include the payoff in the basis! \n");
      exit(-1);
    }

  return(1);
}



// Declaration of allocation/liberation/initialization routines. 
//is made in lmm_bermudaprice.h. Definitions are at the end of this file*/


int computeBermudeanSwaption(long numberMCPaths,int numberTimeStep,
				 Libor *ptLib, Swaption* ptSwpt, Volatility* ptVol, 
				 int RegrBasisDimension,
			          double payoff_as_regressor,
				 char* Basis_Choice,char* Measure_Choice,
			          char Explanatory){
  
  int i,j,l;
  long k;
  double time,dt,AuxNumSpot,AuxOption,AuxScal,p=0.0;
  double *X=NULL,*Res=NULL,*FP=NULL; //X is the state vector for regression
  double *Brownian=NULL,*SwapPrices=NULL,*Numeraire=NULL;
  double *W=NULL;
  int s,numberOfExerciseDates,PayOff_As_Regressor;
  int RegrVarDimension;
  char auxstring[10];
  char ErrorMessage[100];
  Libor *ptLibTemp=NULL;
  Libor *ptLibOld;
  RandomGenerator *ptRand=NULL;
  strcpy(ErrorMessage,"");
  
  //Initialization of Auxiliary Constants:time step for SDE discretization 
  //(Euler scheme), index of first exercing date, payoffasregressor.....
  dt=ptLib->tenor/(double)numberTimeStep;
  s=(int)(ptSwpt->swaptionMaturity/ptLib->tenor);
  numberOfExerciseDates=ptLib->numberOfMaturities-s;
  PayOff_As_Regressor=(int)(payoff_as_regressor/ptLib->tenor)-s;
  
  if ((PayOff_As_Regressor<0)|(PayOff_As_Regressor>(numberOfExerciseDates-1))){
    strcat(ErrorMessage,"Typed payoff_as_regressor value is not valid\n");}
  else if(PayOff_As_Regressor==(numberOfExerciseDates-1)) strcat(ErrorMessage,"Payoff is Not included in regression");
  if((int)(ptSwpt->swapMaturity/ptLib->tenor)!=ptLib->numberOfMaturities){
    strcat(ErrorMessage,"Mismatch between swap maturity and numberof maturities");}

  //RegrVarDim,nametobasis,nametomeasure
  if(Explanatory=='B') RegrVarDimension=ptVol->numberOfFactors;
  else RegrVarDimension=1;
  sprintf(auxstring,"%d",RegrVarDimension);
  strcat(Basis_Choice,auxstring);//Basis_Choice Must be "CanD$(Regr_Var_Dimension)" or "HerD$(Regr_Var_Dimension)"
  Name_To_Basis(ErrorMessage, Basis_Choice,&Basis,RegrVarDimension);
  Name_To_Measure(ErrorMessage, Measure_Choice,&computeEvolution);

  //Bermuda Variables Memory Allocation
  mallocLibor(&ptLibOld,ptLib->numberOfMaturities,ptLib->tenor);//ptLibOld keeps record of the initial values
  mallocLibor(&ptLibTemp,ptLib->numberOfMaturities,ptLib->tenor);
  mallocRandomGenerator(&ptRand,ptVol->numberOfFactors);  
  mallocBermudaVar(&X,&Res,&Brownian,&SwapPrices,&Numeraire,&FP,numberMCPaths,RegrVarDimension,RegrBasisDimension,numberOfExerciseDates,ptVol->numberOfFactors,ErrorMessage);
  W=(double*)malloc(ptVol->numberOfFactors*sizeof(double));
  
  
  
  //Libor Time evolution + record of Brownian Paths, Numeraire Paths and Swap Prices 
  for(k=0;k<numberMCPaths;k++){
    time=0.0;
    AuxNumSpot=(1.0+ptLib->tenor*ptLib->libor[0]);//AuXNumSpot=NumerSpot(T_1)
    copyLibor(ptLibOld,ptLib);
    Set_to_Zero(W,ptVol->numberOfFactors);
  
    for (j=0;j<=(ptLib->numberOfMaturities-2);j++){
      for(i=0;i<numberTimeStep;i++){
	//time evolution from T_j+i*dt to T_j+(i+1)*dt 
	randomVector(ptRand);
	computeEvolution(ptRand,ptLib,ptLibTemp,ptVol,dt,time);
	copyLibor(ptLibTemp,ptLib); 
	time+=dt;
	for(l=0;l<ptVol->numberOfFactors;l++) W[l]+=(sqrt(dt)*ptRand->val[l]);
      }
      AuxNumSpot*=(1.0+ptLib->tenor*ptLib->libor[j+1]);//AuxNumSpot=NumeraireSpot(T_{j+2})
      computeNumeraire(Measure_Choice,ptLib,ptSwpt,Numeraire,j,k,AuxNumSpot);
      if ((s-1)<=j) {
	SwapPrices[k*numberOfExerciseDates+(j+1-s)]=computeSwapPrice(ptLib,ptSwpt,j+1,j+1,ptLib->numberOfMaturities);
	for(l=0;l<ptVol->numberOfFactors;l++) Brownian[k*(numberOfExerciseDates*ptVol->numberOfFactors)+(j+1-s)*ptVol->numberOfFactors+l]=W[l];
	
      }
      putLiborToZero(ptLib,j+1);
      
    }
    
  }
   
  //Backward programming  
  
  //Price at last Excercise date is just excercise price  
  for(k=0;k<numberMCPaths;k++){
    FP[k]=ppos(SwapPrices[k*numberOfExerciseDates+(numberOfExerciseDates-1)]);
    //If s=e-1 the option is indeed European, actualization changes....
    if(numberOfExerciseDates==1) {FP[k]/=Numeraire[k*numberOfExerciseDates+0];}
    else {FP[k]*=(Numeraire[k*numberOfExerciseDates+(numberOfExerciseDates-2)]/Numeraire[k*numberOfExerciseDates+(numberOfExerciseDates-1)]);}
  }
  //Price at time T_j: regression of FP over X(T_j))
  for(j=numberOfExerciseDates-2;j>=0;j--){
    initStateVector(X,Brownian,SwapPrices,Numeraire,j,numberMCPaths,numberOfExerciseDates,ptVol->numberOfFactors,Explanatory);
    Regression(numberMCPaths,numberOfExerciseDates,RegrBasisDimension,ptVol->numberOfFactors,j,
	       PayOff_As_Regressor,X,FP,SwapPrices,Res,ErrorMessage );
    
    for(k=0;k<numberMCPaths;k++){//exercise value
      AuxOption=ppos(SwapPrices[k*numberOfExerciseDates+j]);
      //approximated continuation value, only the at-the-monney paths are taken into account
      if (AuxOption>0){
	// if PayOff_As_Regressor<=j, excercise value is introduced into  regression basis
	if (PayOff_As_Regressor<=j){
	  AuxScal=AuxOption*Res[0];
	  for (l=1;l<RegrBasisDimension;l++){
	    AuxScal+=Basis(X+k*RegrVarDimension,l-1)*Res[l];
	  }
	} else {
	  AuxScal=0.;
	  for (l=0;l<RegrBasisDimension;l++){
	    AuxScal+=Basis(X+k*RegrVarDimension,l)*Res[l];
	  }
	}
	// AuxScal contains the approximated continuation value
	// if AuxScal< exercise value, the optimal stopping time is modified
	if (AuxOption>AuxScal){
	  FP[k]=AuxOption;
	}
      }
      //Discount Factor from time T_{s+j} to T_{s+j-1}}
      if(j>0) FP[k]*=(Numeraire[k*numberOfExerciseDates+j-1]/Numeraire[k*numberOfExerciseDates+j]);
      //Discount from T_s downto  time=0.0 
      else FP[k]*=(1.0/Numeraire[k*numberOfExerciseDates+0]);
    }

  }

   
  //price at time t=0.0 is only the mean discounted value
  ptSwpt->price=0.0;
  for(k=0;k<numberMCPaths;k++){
    ptSwpt->price+=FP[k];
  }
  ptSwpt->price/=(double)numberMCPaths;

  
  if(ErrorMessage!="") printf("\n WARNING: following errors found\n%s\n",ErrorMessage);
  //freeing memory
  freeBermudaVar(&X,&Res,&Brownian,&SwapPrices,&Numeraire);
  freeLibor(&ptLibTemp);
  freeLibor(&ptLibOld);
  freeRandomGenerator(&ptRand);
  free(W);
  return(1);
}

/***************************************************************/
/*                  Auxiliary routines                         */
/***************************************************************/

void mallocBermudaVar(double **RegrVar,double **Res,double** Brownian,double** SwapPrices,double** Numeraire,double **FP,
			long numberMCPaths, int RegrVarDimension,int RegrBasisDimension,int numberOfExerciseDates,int Brown_factors,
			char* ErrorMessage){
  if ((*RegrVar)==NULL){
    (*RegrVar)=(double*)malloc(numberMCPaths*RegrVarDimension*sizeof(double));
    if ((*RegrVar)==NULL) strcat(ErrorMessage,"MEMORY ALLOCATION ERROR,bermudaprice.c\n" );
  }
  if ((*Res)==NULL){
    (*Res)=(double*)malloc(RegrBasisDimension*sizeof(double));
    if ((*Res)==NULL) strcat(ErrorMessage,"MEMORY ALLOCATION ERROR,bermudaprice.c\n" );
  }
  if ((*Brownian)==NULL){
    (*Brownian)=(double*)malloc(numberMCPaths*numberOfExerciseDates*Brown_factors*sizeof(double));
    if ((*Brownian)==NULL) strcat(ErrorMessage,"MEMORY ALLOCATION ERROR,bermudaprice.c\n" );
  }
  if ((*SwapPrices)==NULL){
    (*SwapPrices)=(double*)malloc(numberMCPaths*numberOfExerciseDates*sizeof(double));
    if ((*SwapPrices)==NULL) strcat(ErrorMessage,"MEMORY ALLOCATION ERROR,bermudaprice.c\n" );
  }
  if ((*Numeraire)==NULL){
    (*Numeraire)=(double*)malloc(numberMCPaths*numberOfExerciseDates*sizeof(double));
    if ((*Numeraire)==NULL) strcat(ErrorMessage,"MEMORY ALLOCATION ERROR,bermudaprice.c\n" );
  }
  if ((*FP)==NULL){
    (*FP)=(double*)malloc(numberMCPaths*sizeof(double));
    if ((*FP)==NULL) strcat(ErrorMessage,"MEMORY ALLOCATION ERROR,bermudaprice.c\n" );
  }
 }
void freeBermudaVar(double **RegrVar,double **Res,double** Brownian,double** SwapPrices,double** Numeraire){
  
  free(*RegrVar);
  free(*Res);
  free(*Brownian);
  free(*SwapPrices);
  free(*Numeraire);
  (*RegrVar)=NULL;
  (*Res)=NULL;
  (*Brownian)=NULL;
  (*SwapPrices)=NULL;
  (*Numeraire)=NULL;
  return;
}
/*initialization of the Explicatory variable: Paying Value, Numeraire or Brownian Motion*/ 
void initStateVector(double *X,double *Brownian,double* SwapPrices,double *Numeraire,int j,long numberMCpaths,int numberOfExerciseDates,int numberOfFactors,char Rflag){
  int k,l;
  
  for(k=0;k<numberMCpaths;k++){
    if (Rflag=='B'){for(l=0;l<numberOfFactors;l++) X[k*numberOfFactors+l]=Brownian[k*(numberOfFactors*numberOfExerciseDates)+j*numberOfFactors+l];
    }
    else if(Rflag=='S')  {X[k]=ppos(SwapPrices[k*numberOfExerciseDates+j]);}
    else  X[k]=Numeraire[k*numberOfExerciseDates+j];
			   
  } 
  return;
}    
//******************************************************************************/
//Regression FOR Bermudan Swaptions.
//It is of type  
//     FP~scalprod(Res,VBase(X))
//with FP's,X's as inputs, Res as output. 
//Least squares approach: 
//Res=argmin_a[sum_k (FP^k-scalprod(a,VBase(X^k)))^2]
//Solution obtained by differentiating wrt to Res, rewriting as M*Res=AuxR and then 
//finding Cholesky square root of M.
//We regress only on at the money path (excercing value are needed by algorithm) 
/*******************************************************************************/
void Regression(long NumberMCPaths,int numberOfExerciseDates,
		       int RegrBasis_Dimension, 
		       int X_Dimension, 
		       int Swap_Entry_Time,
		       int PayOff_As_Regressor,
		       double* X,double* FP,double* Swap_Prices,
		       double* Res,
		       char* ErrorMessage )
{
  int i,j;
  double *XPaths_ptr=X,*imDimApppj,AuxOption,AuxScal=0.0;
  double *M=NULL,*AuxR=NULL,*VBase=NULL;
  long k,InTheMoney=0;
  char AuxString[100];


  //Memory Allocation for auxiliary  pointers
  if (M==NULL){
    M=(double*)malloc(RegrBasis_Dimension*RegrBasis_Dimension*sizeof(double));
    if (M==NULL) strcat(ErrorMessage,"MEMORY ALLOCATION ERROR, M, bermudaprice.c\n" );
  }
  if (AuxR==NULL){
    AuxR=(double*)malloc(RegrBasis_Dimension*sizeof(double));
    if (AuxR==NULL) strcat(ErrorMessage,"MEMORY ALLOCATION ERROR, AuxR, bermudaprice.c\n" );
  }
   if (VBase==NULL){
    VBase=(double*)malloc(RegrBasis_Dimension*sizeof(double));
    if (VBase==NULL) strcat(ErrorMessage,"MEMORY ALLOCATION ERROR, VBase, bermudaprice.c\n" );
  }
  //Initialisation of auxiliary pointers
  for (i=0;i<RegrBasis_Dimension;i++){
	AuxR[i]=0;
    for (j=0;j<RegrBasis_Dimension;j++) M[i*RegrBasis_Dimension+j]=0;
  }
  for(k=0;k<NumberMCPaths;k++){
  //kth regressor value 
	AuxOption=ppos(Swap_Prices[k*numberOfExerciseDates+Swap_Entry_Time]);
	//only the at-the-monney path are taken into account
	if (AuxOption>0){
	  InTheMoney++;
	  //value of the regressor basis on the kth path
	  if (PayOff_As_Regressor<=Swap_Entry_Time){
	    //here, the payoff function is introduced in the regression basis
		VBase[0]=AuxOption;
		for (i=1;i<RegrBasis_Dimension;i++){
		  VBase[i]=Basis(XPaths_ptr,i-1);
		}
	  } else {
		for (i=0;i<RegrBasis_Dimension;i++){
		  VBase[i]=Basis(XPaths_ptr,i);
		}
	  }
	  //empirical regressor dispersion matrix
	  imDimApppj=M;
	  for (i=0;i<RegrBasis_Dimension;i++)
		for (j=0;j<RegrBasis_Dimension;j++){
		  *imDimApppj+=VBase[i]*VBase[j];
		  imDimApppj++;
		}
	  //auxiliary for regression formulae
	  for (i=0;i<RegrBasis_Dimension;i++){
		AuxR[i]+=FP[k]*VBase[i];  
	  }
    }
	XPaths_ptr+=X_Dimension;
  }
  if (InTheMoney==0){
    strcat(ErrorMessage,"Warning: there are NO InTheMoney Paths");
    sprintf(AuxString,"%d",Swap_Entry_Time);
    strcat(ErrorMessage,AuxString);
    for (i=0;i<RegrBasis_Dimension;i++){Res[i]=0;}
  } 
  else {
    if (InTheMoney<NumberMCPaths/3){
      strcat(ErrorMessage,"Warning: there are very few InTheMoney Paths");
      sprintf(AuxString,"%d (%ld out of %ld)",Swap_Entry_Time,InTheMoney,NumberMCPaths);
      strcat(ErrorMessage,AuxString);
    }
    //normalisation
    for (i=0;i<RegrBasis_Dimension;i++){
      AuxR[i]/=(double)InTheMoney;
      for (j=0;j<RegrBasis_Dimension;j++){M[i*RegrBasis_Dimension+j]/=(double)InTheMoney;}
    }
    //Cholesky algorithm
    if (Cholesky(M,RegrBasis_Dimension)){
      strcat(ErrorMessage," Cholesky_Warning ") ;
      sprintf(AuxString,"%d",Swap_Entry_Time);
      strcat(ErrorMessage,AuxString);
    }
    //computation of the regression coefficients
    Resolution(AuxR,Res,M,RegrBasis_Dimension);
  }

  free(M);
  free(AuxR);
  free(VBase);
}
