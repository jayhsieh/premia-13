//Xiao WEI(CUFE Bejing) and Nicolas PRIVAULT(University of Poitiers) 2006
//Calibrate eta1, eta2, rhoinf, b and ginf jointly of Libor model to swaptions
//Schoenmakers Algorithm : "Calibration of LIBOR models to caps and swaptions: a way arounds 
//intrinsic instabilities via parsimonious structures and collateral criterion".

//Version modified by Xiao in 30 August, 2008

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include "min.h"
#include "f2c.h"
#include <cstring>
#include <string>
#include <ctype.h>
#include <cctype>
using namespace std;
using std::string;
#define MAXLINE 100
#define DELTA 0.5

double * B, * Bi, * Gamma, * W, * Y, * L, * T, * Shift, * ShiftMark, ** SigmaB, ** Smile;// pointers for market data memory
int * Mat, * Exp; // pointers for maturities in tenor structure
double l[5], u[5], x[5], ParaValue[5];// lower bounds, upper bounds, initial values, renewed values of calibrate parameters by user
int NUM_VALS_P, MAX_TENOR, NUM_SWAP;//up to maturity of market data used to calibrate, maximum number of tenor structure in given market data, number of given swaption volatilities
int Maxmatu, Maxexp, MaxMark, MiddleMark;// maximum number of maturity, maximum number of expirity, maximum number of shift in smile, middle number of shift in smile
fpos_t pos;
double delta = 0.5;


void ChooseCaliParameter(char ParaName[100], double LeftBound, double RightBound, double InitialValue, int i)
{
	string a;
	cout<<endl<<"________Parameter  "<<ParaName <<"?"<<endl;
	cout<<"		"<<"(Calibrate: y, No: n)   ? ";
	cin.clear();
	getline(cin,a, '\n');
	if (a=="n") 
	{
		loop2:
		cout<<"Not calibrate parameter "<<ParaName<<endl;
		cout<<ParaName<<"= "<<InitialValue<<endl;
		cout<<"		"<<"(Ok: y, Modify: m)   ? ";
		cin.clear();
		getline(cin,a, '\n');
		if (a=="m")
		{
			cout<<"		New value:          ? ";
			cin >> ParaValue[i];
			cin.ignore();
			if (ParaValue[i] <= RightBound && ParaValue[i] >= LeftBound)
			{
				u[i] = ParaValue[i];
				l[i] = ParaValue[i];
				x[i] = ParaValue[i];
			}
 			else
			{
				printf("\n		Error: Invalid value! Input %s satisfied %f<=%s<=%f !\n",ParaName,LeftBound,ParaName,RightBound);
				printf("\n");
				goto loop2;
			}
		}
		else 
		{
			u[i]=InitialValue;
			l[i]=InitialValue;
			x[i]=InitialValue;
		}
	}
	else 
	{
		loop1:
		cout<<"Calibrate parameter "<<ParaName<<endl;
		cout<<"Initial value: "<<InitialValue<<endl;
		cout<<"Lower bound:   "<<LeftBound<<endl;
		cout<<"Upper bound:   "<<RightBound<<endl;

		
		cout<<"		"<<"(Ok: y, Modify: m)    ? ";
		cin.clear();
		getline(cin,a, '\n');
		if (a=="m")
		{
			double testupper, testlower, testinitial;	
			cout<<"		Initial value:   ? ";
			cin>>testinitial;
			cin.ignore();
			cout<<"		Lower bound:     ? ";
			cin>>testlower;
			cin.ignore();
			cout<<"		Upper bound:     ? ";
			cin>>testupper;
			cin.ignore();
			if (RightBound>=testupper && testupper>=testinitial && testinitial>=testlower && testlower>=LeftBound)
			{
				u[i]= testupper;
				l[i]= testlower;
				x[i]= testinitial;
			}
			else
			{
				cout<<"			Error: invalid value! "<<endl;
				goto loop1;
			}
		}
		else 
		{
			u[i] = RightBound;
			l[i] = LeftBound;
			x[i] = InitialValue;
		}
	}
}
void ReadMark(FILE *fpt,int marksize,char markname[100])
{
	int errormark;
	signed char ch;
	char tchar[100];
	int nn=0;
	int mm=0;
	if(fpt == NULL)
	{
		printf("\n Error: File couldn't be opened (1).\n");
		exit(-1);
	}
	else
	{
  		if (errormark!=1)
		{
			for (int i=1; i<1000000;i++) 
 			{
				while(ch!='<')
				{
 					ch = getc(fpt);
					if (ch==EOF)
					{
						cout<<"Error: Couldn't find  " << markname<<">  in the XML file"<<endl;
						exit(1);
					}
 				}
 				if (ch=='<')
				{
					do 
					{
						tchar[nn]=ch;
						ch = getc(fpt);
						if (tchar[nn]==markname[nn])
						{
							nn++;
							if (nn == (marksize-1))
							{
								break;
							}
						}
						else
						{
							tchar[nn]=ch;
							nn=0;
						}
				
					}while(ch!='>');
				}
 				if (nn == (marksize-1) )
				{
					break;
				}
				else
				{
					ch = getc(fpt);
					if (ch==EOF)
					{
						printf("\nThere is something wrong to find %s>\n",markname);
						errormark=1;
						exit(1);
					}
				}
			}
			fgetpos (fpt,&pos);
		}
	}
}
void ReadMaturity(FILE *fpt)
{
	char teststr[30];
	fsetpos (fpt,&pos);
	int mm=0;
	do 
	{
		fscanf(fpt, "%s", &teststr);
		Mat[mm]=atoi(teststr);
		mm++;
	}while ((atoi(teststr))!=0);
	Maxmatu = mm-1;
}
void ReadExpiry(FILE *fpt)
{
	char teststr[30];
	fsetpos (fpt,&pos);
	int mm=0;
	do 
	{
		fscanf(fpt, "%s", &teststr);
		Exp[mm]=atoi(teststr);
		mm++;
	}while ((atoi(teststr))!=0);
	Maxexp = mm-1;
}
void ReadSwapRate(FILE *fpt,double vect[100])// read discount factor or forward libor rate 
{
	char teststr[30];
	fsetpos (fpt,&pos);
	int mm=0;
	do 
	{
		fscanf(fpt, "%s", &teststr);
		vect[(Mat[mm])]= atof(teststr);
		mm++;
	}while (mm<Maxmatu);
}
void ReadCapVola(FILE *fpt,double vect[100])
{
	char teststr[30];
	fsetpos (fpt,&pos);
	int i=0;
	do 
	{
		fscanf(fpt, "%s", &teststr);
		Gamma[(Mat[i])]=0.01*atof(teststr);
		i++;
	}while (i<Maxmatu);
}
void ReadShift(FILE *fpt)
{
	char teststr[30];
	fsetpos (fpt,&pos);
	int i=0;
	do 
	{
		fscanf(fpt, "%s", &teststr);
		Shift[(Mat[i])]=0.01*atof(teststr);
		i++;
	}while (i<Maxmatu);
}
void ReadShiftMark(FILE *fpt, char EndMark[100], double vect[100])
{
	char teststr[30];
	fsetpos (fpt,&pos);
 	int i=0;
	do 
	{
		i++;
		fscanf(fpt, "%s", &teststr);
		vect[i]=0.01*atof(teststr);
	}while(atof(teststr)!=0);
	MiddleMark=i;
	do
	{
		i++;
		fscanf(fpt, "%s", &teststr);
		vect[i]=0.01*atof(teststr);
	}while(atof(teststr)!=0);
	MaxMark=i;
}
void ReadSwapVola(FILE *fpt)
{	
	char teststr[30];
	fsetpos (fpt,&pos);
	int mm=0;
	int nn=0;
	do
	{
		do 
		{
			fscanf(fpt, "%s", &teststr);
			SigmaB[Mat[mm]][Mat[mm]+Exp[nn]]=0.01*atof(teststr);
			nn++;
		}while (nn<Maxexp);
		mm++;
		nn=0;
	}while(mm<Maxmatu);
}
void ReadCapSmile(FILE *fpt)
{
	char teststr[30];
	fsetpos (fpt,&pos);
	int mm=0;
	int nn=0;
	do
	{
		for (int i=1; i<MiddleMark; i++) 
		{
			fscanf(fpt, "%s", &teststr);
			Smile[i][Mat[mm]]=(-1)*0.01*atof(teststr);
		}
		for (int i=MiddleMark; i<MaxMark; i++) 
		{
			fscanf(fpt, "%s", &teststr);
			Smile[i][Mat[mm]]=0.01*atof(teststr);
		}
		mm++;
		nn=0;
	}while(mm<Maxmatu);
}
void CapItmOtmToAtm()
{
	double k;
	double Difference;
	int Mark[100];
	for (int j=1; j<MAX_TENOR-1;j++) 
	{
		Mark[j]=1;
		Difference=20.0;
		for (int i=1; i< MaxMark; i++)
		{
			if ((fabs(Shift[j]-ShiftMark[i]))< Difference)	
			{
				Difference = fabs(Shift[j]-ShiftMark[i]);
				Mark[j]=i;
			}
		}
		k = Gamma[j];
		Gamma[j]= k/(1+Smile[(Mark[j])][j]);
	}
}
void SwapItmOtmToAtm()
{
 	double k;
	double Difference;
	int Mark[100];
	for (int j=1; j<MAX_TENOR-1;j++) 
	{
		Mark[j]=1;
		Difference=20.0;
		for (int i=1; i< MaxMark; i++)
		{
			if ((fabs(Shift[j]-ShiftMark[i]))< Difference)	
			{
				Difference = fabs(Shift[j]-ShiftMark[i]);
				Mark[j]=i;
			}
		}
		for (int i=1; i<MAX_TENOR;i++)		
		{
			if (SigmaB[j][j+i]!=0.0)
			{
				k= SigmaB[j][i];
				SigmaB[j][i]=k/(1+Smile[(Mark[j])][j]);
			}
		}
	}
}

void Interpolation(int vecsize,double vect[100])
{
	int j = 1;
	double aver;
	int * mark= new int [MAX_TENOR]; 
	for (int i =1; i<vecsize;i++)
	{
		if(vect[i]!=0.0)
		{
			mark[j]=i;
			j++;
		}
	}
	int maxmark = j;
	mark[0]=1;
	mark[maxmark]=vecsize-1;
	for (int j=0;j<maxmark;j++)
	{
		if (vect[(mark[j])]!=0.0 && vect[(mark[(j+1)])]!=0.0)
		{
			aver = (vect[(mark[j+1])]-vect[(mark[j])])/(mark[j+1]-mark[j]);
			for (int i=mark[j]; i<mark[(j+1)];i++)
			{
				vect[i]=vect[(mark[j])]+(i-mark[j])*aver;
			}
		}
		else if (vect[(mark[j+1])]==0.0)
		{
			aver = (vect[(mark[j])]-vect[(mark[j-1])])/(mark[j]-mark[j-1]);
			for (int i=mark[j]; i<(mark[(j+1)]+1);i++)
			{
				vect[i]=vect[(mark[j])]+(i-mark[j])*aver;
			}
		}
		else if (vect[mark[j]]==0.0)
		{
			aver = (vect[(mark[j+2])]-vect[(mark[j+1])])/(mark[j+2]-mark[j+1]);
			for (int i=mark[j]; i<mark[(j+1)];i++)
			{
				vect[i]=vect[(mark[j+1])]-(mark[j+1]-i)*aver;
			}
		}
	}
}
void ModifiedSmile(double test, double vect[100])
{
	if(test<0)
	{
		for (int j= 0; j< MAX_TENOR; j++)
		{
			if (vect[j]>0)
			{
				vect[j]=0;
			}
		}
	}
	else
	{
		for (int j= 0; j< MAX_TENOR; j++)
		{
			if (vect[j]<0)
			{
				vect[j]=0;
			}
		}
	}
}
void CalculateDisFac(int vecsize)
{
	double Bii[100];
	Bii[0]=1;
	for (int i=1;i<vecsize; i++)
	{
		double sum = 0.0;
 		for (int j=1; j<i;j++)
 		{
 			sum = sum + Bi[j];
 		}
 		Bii[i]=(1-delta*0.01*B[i]*sum)/(1+0.01*B[i]*delta);
		Bi[i]=Bii[i-1]/(1+0.01*B[i]*delta);
  		printf("swaprate[%d]=%f,Bi[%d]=%f,Bii[%d]=%f\n",i,B[i],i,Bi[i],i,Bii[i]);
	}
}
void PrintValue(int vecsize,double vect[100])
{
	for (int i=1; i<vecsize;i++)
	{
		printf("vect[%d]=%f\n",i,vect[i]);
	}
}
void DefineT()
{
	for (int i=1; i< MAX_TENOR; i++)
	{
		T[i]=i*DELTA;
	}
}

double IntegrateGsDenominator(double b, double ginf, int I)
{
	return double(pow(ginf,2)*T[I] + 2.0*ginf*(1.0-ginf)*(1.0/b)*(1.0-exp((-1)* b *T[I]))+pow((1.0-ginf),2)*(1.0/(2.0*b))*(1.0-exp((-2)*b*T[I])));
}
double IntegrateGsDenominatorMSF (double b, double ginf, int i, int p)
{
	return double(pow(ginf,2)*T[p] + 2.0*ginf*(1.0-ginf)*(1.0/b)*(exp(b*(T[p]-T[i]))-exp((-1)*b*T[i])) + pow((1.0-ginf),2)*(1.0/(2.0*b))*(exp(2.0*b*(T[p]-T[i]))-exp((-1)*2.0*b*T[i])));
}
double IntegrateGsNumerator(double b, double ginf, int I, int J, int p)
{
	return double(pow(ginf,2)*T[p]+pow((1.0-ginf),2)*exp((-1)*b*(T[I]+T[J]))*(1.0/(2.0*b))*(exp(2*b*T[p])-1.0)+ginf*(1.0-ginf)*(exp((-1)*b*T[I])+exp((-1)*b*T[J]))*(1.0/b)*(exp(b*T[p])-1.0));
}
double CalculateFraction(double b, double ginf, int I , int J , int p)
{
	return double(IntegrateGsNumerator(b,ginf,I,J,p)/(sqrt(IntegrateGsDenominator(b,ginf,I)*IntegrateGsDenominator(b,ginf,J))));
}
double CalculateAlpha(double b, double ginf, int I , int J , int p)
{
	return double((sqrt(T[I]*T[J])/T[p])*CalculateFraction(b,ginf,I,J,p));
}
double CalculateAlphaMSF(double b, double ginf, int i, int j, int p)
{
	return double (IntegrateGsNumerator(b,ginf,i,j,p)/sqrt(IntegrateGsDenominatorMSF(b, ginf, i, p)*IntegrateGsDenominatorMSF(b, ginf,j, p)));
}
double CalculateC1(int i, int j)
{
	return double ((i*i+j*j+i*j-3.0*(MAX_TENOR-1)*i-3.0*(MAX_TENOR-1)*j+3.0*i+3.0*j+2.0*(MAX_TENOR-1)*(MAX_TENOR-1)-(MAX_TENOR-1)-4.0)/(((MAX_TENOR-1)-2.0)*((MAX_TENOR-1)-3.0)));
}
double CalculateC2(int i, int j)
{
	return double ((i*i+j*j+i*j-(MAX_TENOR-1)*i-(MAX_TENOR-1)*j-3.0*i-3.0*j+3.0*(MAX_TENOR-1)+2.0)/(((MAX_TENOR-1)-2.0)*((MAX_TENOR-1)-3.0)));
}
double CalculateRho(double eta1, double eta2, double rhoinf, int i, int j)
{
	double rho = 0.0;
	if (0.0 <=eta2 && eta2<=3.0*eta1 && (eta1+eta2)<=-log(rhoinf))
	{
		rho =  exp(-fabs(-j+i)/((MAX_TENOR-1)-1.0)*(-log(rhoinf)+eta1*CalculateC1(i,j)-eta2*CalculateC2(i,j))); 
	}
	return rho;
}
double CalculatePartialSum(int p , int q)
{
	double PartialSum = 0.0;
	for(int i = 1 ; i < (q-p)/2+1; i++)
	{
		PartialSum +=  2*DELTA*Bi[p+2*i];
	}
	return PartialSum;
}
void CalculateL()
{
	for(int i = 1; i < MAX_TENOR; i++)
	{	
		L[i] = (1/DELTA)*(Bi[i]-Bi[i+1])/Bi[i+1];
	}
}
double CalculateF(int i, int q)
{
	double F = 0.0;
	CalculateL();
	for (int j=i; j< q; j++)
	{
		F += DELTA*Bi[j+1]*L[j];
	}
	return F;
}
double CalculateG(int s, int q)
{
	double G = 0.0;
	if (q>=s)
	{
		for (int k= 0; k< ((q-s)/2+1); k++)
		{
			G += 2.0*Bi[s+2*k];
		}
	}
	else 
	{
		G = 0.0;
	}
	return G;
}
void CalculateY(int p, int q)
{
	for(int i = p ; i <q ; i++)
	{
		int k = i/2;
		Y[i] = (CalculateF(p,q)*CalculateG(2*k+2,q)-CalculateF(i,q)*CalculateG(p+2,q))/((1.0+DELTA*L[i])*pow(CalculateG(p+2,q),2));
	}
}
void CalculateW(int p, int q)
{
	CalculateY(p,q);
	for(int i = p ; i <q ; i++)
	{
		W[i] = DELTA*Bi[i+1]/CalculatePartialSum(p,q)+Y[i];
	}
}
double CalculateSpq(int p, int q)
{
	return double((Bi[p]-Bi[q])/CalculatePartialSum(p,q));
}
double CalculateSigmaSquare(double b,double ginf,double eta1,double eta2,double rhoinf,int p, int q)
{
	CalculateL();
	CalculateW(p, q);
	double sigma = 0.0;
	for(int i=p; i<q; i++)
	{
		for(int j=p; j<q; j++)
		{
			sigma+= W[i]*W[j]*L[i]*L[j]*Gamma[i]*Gamma[j]*CalculateAlpha(b,ginf,i,j, p)*CalculateRho(eta1,eta2,rhoinf,i,j)/pow(CalculateSpq(p,q),2);
		}
	}
	return sigma;
}
double CalculateSigmaSquareMSF(double b,double ginf,double eta1,double eta2,double rhoinf,int p,int q)
{
	CalculateL();
	CalculateW(p, q);
	double sigma = 0.0;
	for(int i=p; i<q; i++)
	{
		for(int j=p; j<q; j++)
		{
			sigma+= W[i]*W[j]*L[i]*L[j]*Gamma[i]*Gamma[j]*CalculateAlphaMSF(b, ginf, i, j, p)*CalculateRho(eta1,eta2,rhoinf,i,j)/pow(CalculateSpq(p, q),2);
		}
	}
	return sigma;
}	
double CalculateSigma(double b, double ginf, double eta1, double eta2, double rhoinf,int p,int q)
{
	return double (sqrt(CalculateSigmaSquare(b,ginf,eta1,eta2,rhoinf,p,q)));
}
double CalculateSigmaMSF(double b, double ginf,double eta1, double eta2, double rhoinf,int p, int q)
{
	return double (sqrt(CalculateSigmaSquareMSF(b,ginf,eta1,eta2,rhoinf,p,q)));
}
double Calculatediffe(double b, double ginf, double eta1, double eta2, double rhoinf,int p, int q)
{
	return double (pow(((SigmaB[p][q]-CalculateSigma(b,ginf,eta1,eta2,rhoinf,p,q))/SigmaB[p][q]),2));
}
double CalculatediffeMSF(double b, double ginf, double eta1, double eta2, double rhoinf,int p, int q)
{
	return double (pow(((SigmaB[p][q]-CalculateSigmaMSF(b,ginf,eta1,eta2,rhoinf,p,q))/SigmaB[p][q]),2));
}
double CalculateTotalDifferenceInRms(double b,double ginf,double eta1, double eta2, double rhoinf)
{
	double todiffe = 0.0;
	for(int p = 1; p < NUM_VALS_P ; p++)
	{
		for(int q = p+1 ; q < MAX_TENOR; q++)
		{
			 if (SigmaB[p][q]!= 0.0)
			 {
			 	todiffe = todiffe + Calculatediffe(b, ginf, eta1,eta2,rhoinf, p,q);
			 }	
		}
	}
	return todiffe;
}
double CalculateTotalDifferenceInRmsMSF(double b,double ginf,double eta1,double eta2, double rhoinf )
{
	double todiffe = 0.0;
	for(int p = 1; p < NUM_VALS_P ; p++)
	{
		for(int q = p+1; q < MAX_TENOR; q++)
		{
			 if (SigmaB[p][q]!= 0.0)
			 {
			 	todiffe = todiffe + CalculatediffeMSF(b, ginf,eta1,eta2,rhoinf,p,q);
			 }	
		}
	}
	return todiffe;
}
double CalculateRms(double b,double ginf,double eta1, double eta2, double rhoinf)
{
	double rms = sqrt((8.0/(((NUM_SWAP-2))*((NUM_SWAP)-4)))*CalculateTotalDifferenceInRms(b, ginf,eta1,eta2,rhoinf));
	return rms;
}
double CalculateRmsMSF(double b,double ginf, double eta1, double eta2, double rhoinf)
{
	double rms=sqrt((8.0/(((NUM_SWAP-2))*((NUM_SWAP)-4)))*CalculateTotalDifferenceInRmsMSF(b,ginf,eta1,eta2,rhoinf));
	return rms;
}
double CalculateMS(double b,double ginf,double eta1, double eta2,double rhoinf)
{
	return double(1000.0*(pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),2)*sqrt(pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4))));
}

/////////*****************************************************************///////////////////////////////////
double CalculateSigma2(double eta1, double eta2, double rhoinf, int p, int q)
{
	CalculateL();
	CalculateW(p, q);
	double sigmasqr = 0.0;
	for(int i= p; i<q; i++)
	{
		for(int j=p; j<q; j++)
		{
			sigmasqr+= W[i]*W[j]*L[i]*L[j]*Gamma[i]*Gamma[j]*CalculateRho(eta1, eta2, rhoinf, i, j)/pow(CalculateSpq(p, q),2);
		}
	}
	return sqrt(sigmasqr);
}	
double Calculatediffe2(double eta1, double eta2, double rhoinf,int p, int q)
{	
	return double (pow(((SigmaB[p][q]-CalculateSigma2(eta1,eta2,rhoinf,p,q))/SigmaB[p][q]),2));
}	 
double CalculateTotalDifferenceInRms2(double eta1, double eta2, double rhoinf)
{	
	double todiffe = 0.0;
	for(int p = 1; p < NUM_VALS_P ; p++)
	{
		for(int q = p+1 ; q < MAX_TENOR; q++)
		{
			 if (SigmaB[p][q]!= 0.0)
			 {
			 	todiffe = todiffe + Calculatediffe2(eta1,eta2,rhoinf,p,q);
			 }	
		}
	}
	return todiffe;
}
double CalculateRms2(double eta1, double eta2, double rhoinf)
{
	double rms = sqrt((8.0/(((NUM_SWAP)-2)*((NUM_SWAP)-4)))*CalculateTotalDifferenceInRms2(eta1,eta2,rhoinf));
	return rms;
}


///////////////////****************************************************//////////////////////////////////////

//////////////////.........................derivative....................////////////////////////
double GsDeriB(double b,double ginf, double s)
{
	double GsderiB = (1.0-ginf)*((-1)*s)*exp((-1)*b*s);
	return GsderiB;
}
double GsDeriGinf(double b,double ginf, double s)
{
	double GsderiG = 1-exp((-1)*b*s);
	return GsderiG;
}
double DenominatorDeriB (double b, double g, int i)
{
	return double ( 2.0*T[i]*exp(-2.0*b*T[i])-4.0*T[i]*exp(-2.0*b*T[i])*g+2.0*T[i]*exp(-2.0*b*T[i])*g*g+4.0*T[i]*exp(-b*T[i])*g-4.0*T[i]*exp(-b*T[i])*g*g+2.0*g*g*T[i])/b/2.0
	-(1.0+2.0*g-3.0*g*g-exp(-2.0*b*T[i])+2.0*exp(-2.0*b*T[i])*g-exp(-2.0*b*T[i])*g*g-4.0*exp(-b*T[i])*g+4.0*exp(-b*T[i])*g*g-2.0*g*g*log(exp(-b*T[i])))/(b*b)/2.0;
}
double DenominatorDeriBMSF (double b, double ginf, int i, int p)
{
	return double (2.0*ginf*(1.0-ginf)*(1.0/b)*((T[p]-T[i])*exp(b*(T[p]-T[i]))+T[i]*exp((-1)*b*T[i]))-2.0*ginf*(1.0-ginf)*(1.0/pow(b,2))*(exp(b*(T[p]-T[i]))-exp((-1)*b*T[i]))+pow((1.0-ginf),2)*(1.0/(2.0*b))*(2.0*(T[p]-T[i])*exp(2.0*b*(T[p]-T[i]))+2.0*T[i]*exp((-1)*2.0*b*T[i]))-pow((1.0-ginf),2)*(1.0/(2.0*b*b))*(exp(2.0*b*(T[p]-T[i]))-exp((-1)*2.0*b*T[i])));
}
double DenominatorDeriGinf (double b, double g, int i)
{
	return double ((2.0-6.0*g+2.0*exp(-2.0*b*T[i])-2.0*exp(-2.0*b*T[i])*g-4.0*exp(-b*T[i])+8.0*exp(-b*T[i])*g-4.0*g*log(exp(-b*T[i])))/b/2.0);
}
double DenominatorDeriGinfMSF (double b, double ginf, int i,int p)
{
	return double (2.0*ginf*T[p]+(2.0*(1.0-ginf)-2.0*ginf)*(1.0/b)*(exp(b*(T[p]-T[i]))-exp((-1)*b*T[i]))-2.0*(1.0-ginf)*(1.0/(2.0*b))*(exp(2.0*b*(T[p]-T[i]))-exp((-1)*2.0*b*T[i])));
}
double NumeratorDeriB(double b, double g, int i, int j, int p)
{
	double MapleGenVar1 = 0.0;
	double MapleGenVar2 = 0.0;
	double MapleGenVar3 = 0.0;
	double MapleGenVar4 = 0.0;
	MapleGenVar1 = (-2.0*g*T[i]*exp(b*T[i])+2.0*g*g*T[i]*exp(b*T[i])-2.0*g*T[j]*exp(b*T[j])+2.0
*g*g*T[j]*exp(b*T[j])+2.0*T[p]*g*g*exp(b*T[i]+b*T[j])+2.0*T[p]*g*g*b*(T[i]+T[j])*exp(b*T[i]+b*T[j])+2.0*g*(T[p]+
T[i])*exp(b*T[p]+b*T[i])-2.0*g*g*(T[p]+T[i])*exp(b*T[p]+b*T[i])+2.0*g*(T[p]+T[j])*exp(b*T[p]+b*T[j])+2.0*T[p]*exp(
2.0*b*T[p])-4.0*g*T[p]*exp(2.0*b*T[p])-2.0*g*g*(T[p]+T[j])*exp(b*T[p]+b*T[j])+2.0*g*g*T[p]*exp(2.0*b*T[p])
)/b*exp(-b*T[i]-b*T[j])/2.0;
      MapleGenVar3 = -(-g*g-2.0*g*exp(b*T[i])+2.0*g*g*exp(b*T[i])-2.0*g*exp(b*T[j])-1.0+
2.0*g+2.0*g*g*exp(b*T[j])+2.0*T[p]*g*g*b*exp(b*T[i]+b*T[j])+2.0*g*exp(b*T[p]+b*T[i])-2.0*g*g*exp(
b*T[p]+b*T[i])+2.0*g*exp(b*T[p]+b*T[j])+exp(2.0*b*T[p])-2.0*g*exp(2.0*b*T[p])-2.0*g*g*exp(b*T[p]+b*T[j]
)+g*g*exp(2.0*b*T[p]))/(b*b)*exp(-b*T[i]-b*T[j])/2.0;
      MapleGenVar4 = (-g*g-2.0*g*exp(b*T[i])+2.0*g*g*exp(b*T[i])-2.0*g*exp(b*T[j])-1.0+
2.0*g+2.0*g*g*exp(b*T[j])+2.0*T[p]*g*g*b*exp(b*T[i]+b*T[j])+2.0*g*exp(b*T[p]+b*T[i])-2.0*g*g*exp(
b*T[p]+b*T[i])+2.0*g*exp(b*T[p]+b*T[j])+exp(2.0*b*T[p])-2.0*g*exp(2.0*b*T[p])-2.0*g*g*exp(b*T[p]+b*T[j]
)+g*g*exp(2.0*b*T[p]))/b*(-T[i]-T[j])*exp(-b*T[i]-b*T[j])/2.0;
      MapleGenVar2 = MapleGenVar3+MapleGenVar4;
      return double(MapleGenVar1+MapleGenVar2);
}
double NumeratorDeriGinf(double b, double g, int i, int j, int p)
{
	return double((-2.0*g-2.0*exp(b*T[i])+4.0*g*exp(b*T[i])-2.0*exp(b*T[j])+2.0+4.0*g*exp(b*T[j])+
4.0*T[p]*g*b*exp(b*T[i]+b*T[j])+2.0*exp(b*T[p]+b*T[i])-4.0*g*exp(b*T[p]+b*T[i])+2.0*exp(b*T[p]+b*T[j])-2.0
*exp(2.0*b*T[p])-4.0*g*exp(b*T[p]+b*T[j])+2.0*g*exp(2.0*b*T[p]))/b*exp(-b*T[i]-b*T[j])/2.0);
}
double AlphaDeriB(double b, double g, int i, int j, int p)
{
	
	double MapleGenVar1 = 0.0;
	double MapleGenVar2 = 0.0;
	double MapleGenVar3 = 0.0;
	double MapleGenVar4 = 0.0;
	double MapleGenVar5 = 0.0;
	double MapleGenVar6 = 0.0;
	double MapleGenVar7 = 0.0;
	double MapleGenVar8 = 0.0;
	double MapleGenVar9 = 0.0;
	double MapleGenVar10 = 0.0;
	double MapleGenVar11 = 0.0;
	double MapleGenVar12 = 0.0;
	double MapleGenVar13 = 0.0;
	
	MapleGenVar3 = (-2.0*g*T[i]*exp(b*T[i])+2.0*g*g*T[i]*exp(b*T[i])-2.0*g*T[j]*exp(b*T[j])+2.0
*g*g*T[j]*exp(b*T[j])+2.0*T[p]*g*g*exp(b*T[i]+b*T[j])+2.0*T[p]*g*g*b*(T[i]+T[j])*exp(b*T[i]+b*T[j])+2.0*g*(T[p]+
T[i])*exp(b*T[p]+b*T[i])-2.0*g*g*(T[p]+T[i])*exp(b*T[p]+b*T[i])+2.0*g*(T[p]+T[j])*exp(b*T[p]+b*T[j])+2.0*T[p]*exp(
2.0*b*T[p])-4.0*g*T[p]*exp(2.0*b*T[p])-2.0*g*g*(T[p]+T[j])*exp(b*T[p]+b*T[j])+2.0*g*g*T[p]*exp(2.0*b*T[p])
)*exp(-b*T[i]-b*T[j]);
      MapleGenVar4 = 1/b/sqrt((1.0+2.0*g-3.0*g*g-exp(-2.0*b*T[i])+2.0*exp(-2.0*b*T[i]
)*g-exp(-2.0*b*T[i])*g*g-4.0*g*exp(-b*T[i])+4.0*g*g*exp(-b*T[i])-2.0*g*g*log(exp(-b*T[i])))
/(b*b)*(1.0+2.0*g-3.0*g*g-exp(-2.0*b*T[j])+2.0*exp(-2.0*b*T[j])*g-exp(-2.0*b*T[j])*g*g
-4.0*g*exp(-b*T[j])+4.0*g*g*exp(-b*T[j])-2.0*g*g*log(exp(-b*T[j]))));
      MapleGenVar2 = MapleGenVar3*MapleGenVar4;
      MapleGenVar3 = (-g*g-2.0*g*exp(b*T[i])+2.0*g*g*exp(b*T[i])-2.0*g*exp(b*T[j])-1.0+
2.0*g+2.0*g*g*exp(b*T[j])+2.0*T[p]*g*g*b*exp(b*T[i]+b*T[j])+2.0*g*exp(b*T[p]+b*T[i])-2.0*g*g*exp(
b*T[p]+b*T[i])+2.0*g*exp(b*T[p]+b*T[j])+exp(2.0*b*T[p])-2.0*g*exp(2.0*b*T[p])-2.0*g*g*exp(b*T[p]+b*T[j]
)+g*g*exp(2.0*b*T[p]))*(-T[i]-T[j])*exp(-b*T[i]-b*T[j])/b/sqrt((1.0+2.0*g-3.0*g*g-exp(-2.0*b*T[i]
)+2.0*exp(-2.0*b*T[i])*g-exp(-2.0*b*T[i])*g*g-4.0*g*exp(-b*T[i])+4.0*g*g*exp(-b*T[i])-2.0*g
*g*log(exp(-b*T[i])))/(b*b)*(1.0+2.0*g-3.0*g*g-exp(-2.0*b*T[j])+2.0*exp(-2.0*b*T[j])*g-
exp(-2.0*b*T[j])*g*g-4.0*g*exp(-b*T[j])+4.0*g*g*exp(-b*T[j])-2.0*g*g*log(exp(-b*T[j]))));
      MapleGenVar1 = MapleGenVar2+MapleGenVar3;
      MapleGenVar2 = MapleGenVar1;
      MapleGenVar4 = -(-g*g-2.0*g*exp(b*T[i])+2.0*g*g*exp(b*T[i])-2.0*g*exp(b*T[j])-1.0+
2.0*g+2.0*g*g*exp(b*T[j])+2.0*T[p]*g*g*b*exp(b*T[i]+b*T[j])+2.0*g*exp(b*T[p]+b*T[i])-2.0*g*g*exp(
b*T[p]+b*T[i])+2.0*g*exp(b*T[p]+b*T[j])+exp(2.0*b*T[p])-2.0*g*exp(2.0*b*T[p])-2.0*g*g*exp(b*T[p]+b*T[j]
)+g*g*exp(2.0*b*T[p]))*exp(-b*T[i]-b*T[j])/(b*b)/sqrt((1.0+2.0*g-3.0*g*g-exp(-2.0*b*T[i])+
2.0*exp(-2.0*b*T[i])*g-exp(-2.0*b*T[i])*g*g-4.0*g*exp(-b*T[i])+4.0*g*g*exp(-b*T[i])-2.0*g*g
*log(exp(-b*T[i])))/(b*b)*(1.0+2.0*g-3.0*g*g-exp(-2.0*b*T[j])+2.0*exp(-2.0*b*T[j])*g-exp
(-2.0*b*T[j])*g*g-4.0*g*exp(-b*T[j])+4.0*g*g*exp(-b*T[j])-2.0*g*g*log(exp(-b*T[j]))));
      MapleGenVar6 = -(-g*g-2.0*g*exp(b*T[i])+2.0*g*g*exp(b*T[i])-2.0*g*exp(b*T[j])-1.0+
2.0*g+2.0*g*g*exp(b*T[j])+2.0*T[p]*g*g*b*exp(b*T[i]+b*T[j])+2.0*g*exp(b*T[p]+b*T[i])-2.0*g*g*exp(
b*T[p]+b*T[i])+2.0*g*exp(b*T[p]+b*T[j])+exp(2.0*b*T[p])-2.0*g*exp(2.0*b*T[p])-2.0*g*g*exp(b*T[p]+b*T[j]
)+g*g*exp(2.0*b*T[p]))*exp(-b*T[i]-b*T[j])/2.0;
      MapleGenVar8 = 1/b;
      MapleGenVar10 = 1/sqrt(pow(1.0+2.0*g-3.0*g*g-exp(-2.0*b*T[i])+2.0*exp(-2.0*b
*T[i])*g-exp(-2.0*b*T[i])*g*g-4.0*g*exp(-b*T[i])+4.0*g*g*exp(-b*T[i])-2.0*g*g*log(exp(-b*T[i])
),3.0)/(b*b*b*b*b*b)*pow(1.0+2.0*g-3.0*g*g-exp(-2.0*b*T[j])+2.0*exp(-2.0*b*T[j])*g-
exp(-2.0*b*T[j])*g*g-4.0*g*exp(-b*T[j])+4.0*g*g*exp(-b*T[j])-2.0*g*g*log(exp(-b*T[j])),3.0)
);
      MapleGenVar12 = (2.0*T[i]*exp(-2.0*b*T[i])-4.0*T[i]*exp(-2.0*b*T[i])*g+2.0*T[i]*exp(-2.0
*b*T[i])*g*g+4.0*g*T[i]*exp(-b*T[i])-4.0*g*g*T[i]*exp(-b*T[i])+2.0*g*g*T[i])/(b*b)*(1.0+2.0*g-3.0
*g*g-exp(-2.0*b*T[j])+2.0*exp(-2.0*b*T[j])*g-exp(-2.0*b*T[j])*g*g-4.0*g*exp(-b*T[j])+4.0*g*
g*exp(-b*T[j])-2.0*g*g*log(exp(-b*T[j])));
      MapleGenVar13 = -2.0*(1.0+2.0*g-3.0*g*g-exp(-2.0*b*T[i])+2.0*exp(-2.0*b*T[i])*g
-exp(-2.0*b*T[i])*g*g-4.0*g*exp(-b*T[i])+4.0*g*g*exp(-b*T[i])-2.0*g*g*log(exp(-b*T[i])))/(b
*b*b)*(1.0+2.0*g-3.0*g*g-exp(-2.0*b*T[j])+2.0*exp(-2.0*b*T[j])*g-exp(-2.0*b*T[j])*g*g
-4.0*g*exp(-b*T[j])+4.0*g*g*exp(-b*T[j])-2.0*g*g*log(exp(-b*T[j])))+(1.0+2.0*g-3.0*g*g-
exp(-2.0*b*T[i])+2.0*exp(-2.0*b*T[i])*g-exp(-2.0*b*T[i])*g*g-4.0*g*exp(-b*T[i])+4.0*g*g*exp
(-b*T[i])-2.0*g*g*log(exp(-b*T[i])))/(b*b)*(2.0*T[j]*exp(-2.0*b*T[j])-4.0*T[j]*exp(-2.0*b*T[j])*g
+2.0*T[j]*exp(-2.0*b*T[j])*g*g+4.0*g*T[j]*exp(-b*T[j])-4.0*g*g*T[j]*exp(-b*T[j])+2.0*g*g*T[j]);
      MapleGenVar11 = MapleGenVar12+MapleGenVar13;
      MapleGenVar9 = MapleGenVar10*MapleGenVar11;
      MapleGenVar7 = MapleGenVar8*MapleGenVar9;
      MapleGenVar5 = MapleGenVar6*MapleGenVar7;
      MapleGenVar3 = MapleGenVar4+MapleGenVar5;
      
      return double ((sqrt(T[i]*T[j])/T[p])*(MapleGenVar2+MapleGenVar3));
}
double AlphaDeriBMSF(double b, double g, int i, int j, int p)
{
	return double (NumeratorDeriB(b,g,i,j,p)/sqrt(IntegrateGsDenominatorMSF(b,g,i,p)*IntegrateGsDenominatorMSF(b,g,j,p))-(1.0/2.0)*IntegrateGsNumerator(b,g,i,j,p)/pow(IntegrateGsDenominatorMSF(b,g,i,p)*IntegrateGsDenominatorMSF(b,g,j,p),1.5)*(DenominatorDeriBMSF(b,g,i,p)*IntegrateGsDenominatorMSF(b,g,j,p)+DenominatorDeriBMSF(b,g,j,p)*IntegrateGsDenominatorMSF(b,g,i,p)));
}
double AlphaDeriGinf(double b, double g, int i, int j, int p)
{
	return  double (sqrt(T[i]*T[j])*NumeratorDeriGinf(b,g,i,j,p)/T[p]/sqrt(IntegrateGsDenominator(b,g,i)*IntegrateGsDenominator(b,g,j))-sqrt(T[i]*T[j])*IntegrateGsNumerator(b,g,i,j,p)/T[p]/sqrt(pow(IntegrateGsDenominator(b,g,i),3.0)*pow(IntegrateGsDenominator(b,g,j),3.0))*(DenominatorDeriGinf(b,g,i)*IntegrateGsDenominator(b,g,j)+IntegrateGsDenominator(b,g,i)*DenominatorDeriGinf(b,g,j))/2.0);
}
double AlphaDeriGinfMSF(double b, double g, int i, int j, int p)
{
	return double (NumeratorDeriGinf(b,g,i,j,p)/sqrt(IntegrateGsDenominatorMSF(b,g,i,p)*IntegrateGsDenominatorMSF(b,g,j,p))-(1.0/2.0)*IntegrateGsNumerator(b,g,i,j,p)/pow(IntegrateGsDenominatorMSF(b,g,i,p)*IntegrateGsDenominatorMSF(b,g,j,p),1.5)*(DenominatorDeriGinfMSF(b,g,i,p)*IntegrateGsDenominatorMSF(b,g,j,p)+DenominatorDeriGinfMSF(b,g,j,p)*IntegrateGsDenominatorMSF(b,g,i,p)));
}
double RhoDeriEta1(double eta1, double eta2, double rhoinf, int i, int j)
{
	double deri;
	if (0.0 <=eta2 && eta2<=3.0*eta1 && (eta1+eta2)<=-log(rhoinf))
	{
		 deri=(-fabs(-j+i)/((MAX_TENOR-1)-1.0)*CalculateC1(i,j)*CalculateRho(eta1,eta2,rhoinf,i,j));
	}
	return deri ;
}
double RhoDeriEta2(double eta1, double eta2, double rhoinf, int i, int j)
{
	double deri; 
	if (0.0 <= eta2 && eta2<=3.0*eta1 &&  (eta1+eta2)<= -log(rhoinf))
	{
		 deri=(-fabs(-j+i)/((MAX_TENOR-1)-1.0)*CalculateC2(i,j)*CalculateRho(eta1,eta2,rhoinf,i,j));
	}
	return deri ;
}
double RhoDeriRhoinf(double eta1, double eta2, double rhoinf, int i, int j)
{
	double deri;
	if (0.0<=eta2 && eta2<=3.0*eta1 &&   (eta1+eta2)<=-log(rhoinf))
	{
		
		deri = (fabs(-j+i)/(((MAX_TENOR-1)-1.0)*rhoinf)*CalculateRho(eta1,eta2,rhoinf,i,j));
	}
	return deri;
}
double SigmaSqrtDeriB(double b,double ginf, double eta1, double eta2, double rhoinf,int p,int q)
{
	CalculateL();
	CalculateW(p,q);
	double SigmaSqrtDeriB = 0.0 ;
	for(int i=p;i<q;i++)
	{
		for(int j=p; j<q; j++)
		{
			SigmaSqrtDeriB += W[i]*W[j]*L[i]*L[j]*Gamma[i]*Gamma[j]*CalculateRho(eta1,eta2,rhoinf,i, j)*AlphaDeriB(b,ginf,i,j,p)/pow(CalculateSpq(p,q),2);
		}
	}
	return SigmaSqrtDeriB;
}	
double SigmaSqrtDeriBMSF (double b,double ginf,double eta1,double eta2,double rhoinf,int p,int q)
{
	CalculateL();
	CalculateW(p,q);
	double SigmaSqrtDeriB = 0.0 ;
	for(int i=p;i<q;i++)
	{
		for(int j=p; j<q; j++)
		{
			SigmaSqrtDeriB += W[i]*W[j]*L[i]*L[j]*Gamma[i]*Gamma[j]*CalculateRho(eta1, eta2,rhoinf,i,j)*AlphaDeriBMSF(b,ginf,i,j,p)/pow(CalculateSpq(p,q),2);
		}
	}
	return SigmaSqrtDeriB;
}
double SigmaSqrtDeriGinf (double b,double ginf,double eta1,double eta2,double rhoinf,int p,int q)
{
	CalculateL();
	CalculateW(p,q);
	double SigmaSqrtDeriGinf= 0.0;
	for(int i=p;i<q;i++)
	{
		for(int j=p; j<q; j++)
		{
			SigmaSqrtDeriGinf += W[i]*W[j]*L[i]*L[j]*Gamma[i]*Gamma[j]*CalculateRho(eta1,eta2, rhoinf,i,j)*AlphaDeriGinf(b,ginf,i,j,p)/pow(CalculateSpq(p,q),2);
		}
	}
	return SigmaSqrtDeriGinf;
}
double SigmaSqrtDeriGinfMSF(double b,double ginf,double eta1,double eta2,double rhoinf,int p,int q)
{
	CalculateL();
	CalculateW(p,q);
	double SigmaSqrtDeriGinf= 0.0;
	for(int i=p;i<q;i++)
	{
		for(int j=p; j<q; j++)
		{
			SigmaSqrtDeriGinf += W[i]*W[j]*L[i]*L[j]*Gamma[i]*Gamma[j]*CalculateRho(eta1,eta2, rhoinf,i,j)*AlphaDeriGinfMSF(b,ginf,i,j,p)/pow(CalculateSpq(p,q),2);
		}
	}
	return SigmaSqrtDeriGinf;
}	
double SigmaSqrtDeriEta1(double b,double ginf,double eta1,double eta2,double rhoinf,int p, int q)
{
	CalculateL();
	CalculateW(p,q);
	double SigmaSqrtDeriEta1 = 0.0 ;
	for(int i=p;i<q;i++)
	{
		for(int j=p; j<q; j++)
		{
			SigmaSqrtDeriEta1 += W[i]*W[j]*L[i]*L[j]*Gamma[i]*Gamma[j]*CalculateAlpha(b,ginf,i,j,p)*RhoDeriEta1(eta1,eta2,rhoinf,i, j)/pow(CalculateSpq(p,q),2);
		}
	}
	return SigmaSqrtDeriEta1;
}	
double SigmaSqrtDeriEta1MSF(double b,double ginf,double eta1,double eta2,double rhoinf,int p, int q)
{
	CalculateL();
	CalculateW(p,q);
	double SigmaSqrtDeriEta1 = 0.0 ;
	for(int i=p;i<q;i++)
	{
		for(int j=p; j<q; j++)
		{
			SigmaSqrtDeriEta1 += W[i]*W[j]*L[i]*L[j]*Gamma[i]*Gamma[j]*CalculateAlphaMSF(b,ginf,i,j,p)*RhoDeriEta1(eta1,eta2,rhoinf,i, j)/pow(CalculateSpq(p,q),2);
		}
	}
	return SigmaSqrtDeriEta1;
}	
double SigmaSqrtDeriEta2(double b,double ginf, double eta1, double eta2, double rhoinf,int p, int q)
{
	CalculateL();
	CalculateW(p,q);
	double SigmaSqrtDeriEta2 = 0.0 ;
	for(int i=p;i<q;i++)
	{
		for(int j=p; j<q; j++)
		{
			SigmaSqrtDeriEta2 += W[i]*W[j]*L[i]*L[j]*Gamma[i]*Gamma[j]*CalculateAlpha(b,ginf,i,j,p)*RhoDeriEta2(eta1,eta2,rhoinf,i, j)/pow(CalculateSpq(p,q),2);
		}
	}
	return SigmaSqrtDeriEta2;
}	
double SigmaSqrtDeriEta2MSF(double b,double ginf, double eta1,double eta2,double rhoinf,int p,int q)
{
	CalculateL();
	CalculateW(p,q);
	double SigmaSqrtDeriEta2 = 0.0 ;
	for(int i=p;i<q;i++)
	{
		for(int j=p; j<q; j++)
		{
			SigmaSqrtDeriEta2 += W[i]*W[j]*L[i]*L[j]*Gamma[i]*Gamma[j]*CalculateAlphaMSF(b,ginf,i,j,p)*RhoDeriEta2(eta1,eta2,rhoinf,i, j)/pow(CalculateSpq(p,q),2);
		}
	}
	return SigmaSqrtDeriEta2;
}	
double SigmaSqrtDerirhoinf(double b,double ginf,double eta1,double eta2, double rhoinf,int p, int q)
{
	CalculateL();
	CalculateW(p,q);
	double SigmaSqrtDeriRhoinf = 0.0 ;
	for(int i=p;i<q;i++)
	{
		for(int j=p; j<q; j++)
		{
			SigmaSqrtDeriRhoinf += W[i]*W[j]*L[i]*L[j]*Gamma[i]*Gamma[j]*CalculateAlpha(b,ginf,i,j,p)*RhoDeriRhoinf(eta1,eta2,rhoinf,i, j)/pow(CalculateSpq( p,  q),2);
		}
	}
	return SigmaSqrtDeriRhoinf;
}
double SigmaSqrtDerirhoinfMSF(double b,double ginf,double eta1,double eta2,double rhoinf,int p,int q)
{
	CalculateL();
	CalculateW(p,q);
	double SigmaSqrtDeriRhoinf = 0.0 ;
	for(int i=p;i<q;i++)
	{
		for(int j=p; j<q; j++)
		{
			SigmaSqrtDeriRhoinf += W[i]*W[j]*L[i]*L[j]*Gamma[i]*Gamma[j]*CalculateAlphaMSF(b,ginf,i,j,p)*RhoDeriRhoinf(eta1,eta2,rhoinf,i, j)/pow(CalculateSpq(p,q),2);
		}
	}
	return SigmaSqrtDeriRhoinf;
}
double SigmaDeriB (double b, double ginf, double eta1, double eta2, double rhoinf, int p, int q)
{
	return double (1.0/sqrt(CalculateSigmaSquare(b, ginf, eta1, eta2,rhoinf, p, q))*SigmaSqrtDeriB(b,ginf,eta1,eta2, rhoinf,p, q)/2.0);              
}
double SigmaDeriBMSF (double b, double ginf, double eta1, double eta2, double rhoinf, int p, int q)
{
	return double (1.0/sqrt(CalculateSigmaSquareMSF(b, ginf, eta1,eta2,rhoinf, p, q))*SigmaSqrtDeriBMSF(b,ginf,eta1, eta2,rhoinf,p, q)/2.0);              
}
double SigmaDeriGinf (double b, double ginf, double eta1, double eta2, double rhoinf, int p, int q)
{
	return double (1.0/sqrt(CalculateSigmaSquare(b, ginf, eta1,eta2,rhoinf, p, q))*SigmaSqrtDeriGinf(b,ginf,eta1,eta2, rhoinf,p, q)/2.0);              
}
double SigmaDeriGinfMSF (double b,double ginf, double eta1, double eta2, double rhoinf, int p, int q)
{
	return double (1.0/sqrt(CalculateSigmaSquareMSF(b, ginf, eta1,eta2,rhoinf, p, q))*SigmaSqrtDeriGinfMSF(b,ginf,eta1, eta2,rhoinf,p, q)/2.0);              
}
double SigmaDerieta1(double b, double ginf,double eta1, double eta2, double rhoinf, int p, int q)
{
	return double (1.0/sqrt(CalculateSigmaSquare( b, ginf, eta1,eta2,rhoinf,p,q))*SigmaSqrtDeriEta1( b,ginf,eta1, eta2,rhoinf,p, q)/2.0);              
}
double SigmaDerieta1MSF(double b, double ginf,double eta1, double eta2, double rhoinf, int p, int q)
{
	return double (1.0/sqrt(CalculateSigmaSquareMSF( b, ginf, eta1,eta2,rhoinf,p,q))*SigmaSqrtDeriEta1MSF( b,ginf,eta1,eta2, rhoinf,p, q)/2.0);
}
double SigmaDerieta2(double b, double ginf,double eta1, double eta2, double rhoinf, int p, int q)
{
	return double (1.0/sqrt(CalculateSigmaSquare( b, ginf, eta1,eta2,rhoinf,p,q))*SigmaSqrtDeriEta2( b,ginf,eta1, eta2,rhoinf,p, q)/2.0);              
}
double SigmaDerieta2MSF(double b, double ginf,double eta1, double eta2, double rhoinf, int p, int q)
{
	return double (1.0/sqrt(CalculateSigmaSquareMSF( b, ginf, eta1,eta2,rhoinf,p,q))*SigmaSqrtDeriEta2MSF( b,ginf,eta1,eta2, rhoinf,p, q)/2.0);              
}
double SigmaDerirhoinf(double b, double ginf,double eta1,double eta2,  double rhoinf, int p, int q)
{
	return double (1.0/sqrt(CalculateSigmaSquare( b,  ginf, eta1,eta2,rhoinf,p,q))*SigmaSqrtDerirhoinf( b,ginf,eta1,eta2, rhoinf,p,q)/2.0);              
}
double SigmaDerirhoinfMSF(double b,double ginf,double eta1,double eta2, double rhoinf, int p, int q)
{
	return double (1.0/sqrt(CalculateSigmaSquareMSF( b,  ginf, eta1,eta2,rhoinf,p,q))*SigmaSqrtDerirhoinfMSF( b,ginf,eta1, eta2,rhoinf,p,q)/2.0);              
}
double DifferenceDeriB(double b, double ginf, double eta1, double eta2, double rhoinf,int p, int q)
{
	return double (-2.0*(SigmaB[p][q]-CalculateSigma(b,ginf,eta1,eta2,rhoinf,p,q))/(SigmaB[p][q]*SigmaB[p][q])*SigmaDeriB(b,ginf,eta1,eta2,rhoinf,p,q));
}
double DifferenceDeriBMSF(double b,double ginf, double eta1, double eta2, double rhoinf,int p, int q)
{
	return double (-2.0*(SigmaB[p][q]-CalculateSigmaMSF(b,ginf,eta1,eta2,rhoinf,p,q))/(SigmaB[p][q]*SigmaB[p][q])*SigmaDeriBMSF(b,ginf,eta1,eta2,rhoinf,p,q));
}
double DifferenceDeriGinf(double b,double ginf, double eta1, double eta2, double rhoinf,int p, int q)
{
	return double (-2.0*(SigmaB[p][q]-CalculateSigma(b,ginf,eta1,eta2,rhoinf,p,q))/(SigmaB[p][q]*SigmaB[p][q])*SigmaDeriGinf(b,ginf,eta1,eta2,rhoinf,p,q));
}
double DifferenceDeriGinfMSF(double b,double ginf,double eta1,double eta2,double rhoinf,int p, int q)
{
	return double (-2.0*(SigmaB[p][q]-CalculateSigmaMSF(b,ginf,eta1,eta2,rhoinf,p,q))/(SigmaB[p][q]*SigmaB[p][q])*SigmaDeriGinfMSF(b,ginf,eta1,eta2,rhoinf,p,q));
}
double DifferenceDerieta1(double b,double ginf,double eta1, double eta2, double rhoinf,int p, int q)
{
	return double (-2.0*(SigmaB[p][q]-CalculateSigma(b,ginf,eta1,eta2,rhoinf,p,q))/(SigmaB[p][q]*SigmaB[p][q])*SigmaDerieta1(b,ginf,eta1,eta2,rhoinf,p,q));
}
double DifferenceDerieta1MSF(double b,double ginf,double eta1,double eta2,double rhoinf,int p, int q)
{
	return double (-2.0*(SigmaB[p][q]-CalculateSigmaMSF(b,ginf,eta1,eta2,rhoinf,p,q))/(SigmaB[p][q]*SigmaB[p][q])*SigmaDerieta1MSF(b,ginf,eta1,eta2,rhoinf,p,q));
}
double DifferenceDerieta2(double b,double ginf,double eta1,double eta2,double rhoinf,int p,int q)
{
	return double (-2.0*(SigmaB[p][q]-CalculateSigma(b,ginf,eta1,eta2,rhoinf,p,q))/(SigmaB[p][q]*SigmaB[p][q])*SigmaDerieta2(b,ginf,eta1,eta2,rhoinf,p,q));
}
double DifferenceDerieta2MSF(double b,double ginf,double eta1,double eta2,double rhoinf,int p, int q)
{
	return double (-2.0*(SigmaB[p][q]-CalculateSigmaMSF(b,ginf,eta1,eta2,rhoinf,p,q))/(SigmaB[p][q]*SigmaB[p][q])*SigmaDerieta2MSF(b,ginf,eta1,eta2,rhoinf,p,q));
}
double DifferenceDerirhoinf(double b,double ginf,double eta1,double eta2, double rhoinf,int p, int q)
{
	return double (-2.0*(SigmaB[p][q]-CalculateSigma(b,ginf,eta1,eta2,rhoinf,p,q))/(SigmaB[p][q]*SigmaB[p][q])*SigmaDerirhoinf(b,ginf,eta1,eta2,rhoinf,p,q));
}
double DifferenceDerirhoinfMSF(double b,double ginf,double eta1,double eta2,double rhoinf,int p,int q)
{
	return double (-2.0*(SigmaB[p][q]-CalculateSigmaMSF(b,ginf,eta1,eta2,rhoinf,p,q))/(SigmaB[p][q]*SigmaB[p][q])*SigmaDerirhoinfMSF(b,ginf,eta1,eta2,rhoinf,p,q));
}
double TotalDeriB(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	double TotalDiffDeriB = 0.0;
	for( int p = 1 ; p < NUM_VALS_P ; p++)
	{
			for( int q = p+1 ; q < MAX_TENOR; q++)
			{
				if (SigmaB[p][q]!=0.0)
				{
					TotalDiffDeriB += DifferenceDeriB(b, ginf,eta1,eta2,rhoinf, p, q);
				}
			 }
	 }
	return double((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*TotalDiffDeriB);
}
double TotalDeriBMSF(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	double TotalDiffDeriB = 0.0;
	for( int p = 1 ; p < NUM_VALS_P ; p++)
	{
			for( int q = p+1 ; q < MAX_TENOR; q++)
			{
				if (SigmaB[p][q]!=0.0)
				{
					TotalDiffDeriB += DifferenceDeriBMSF(b, ginf,eta1,eta2,rhoinf, p, q);
				}
			 }
	 }
	return double((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*TotalDiffDeriB);
}
double TotalDeriGinf(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	double TotalDiffDeriGinf = 0.0;
	for( int p = 1 ; p < NUM_VALS_P ; p++)
		{
			for( int q = p+1 ; q < MAX_TENOR; q++)
			{
				if (SigmaB[p][q]!=0.0)
				{
					TotalDiffDeriGinf += DifferenceDeriGinf(b, ginf, eta1,eta2,rhoinf,p, q);
				}
			 }
		}
	return double((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*TotalDiffDeriGinf);
}
double TotalDeriGinfMSF(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	double TotalDiffDeriGinf = 0.0;
	for( int p = 1 ; p < NUM_VALS_P ; p++)
		{
			for( int q = p+1; q < MAX_TENOR; q++)
			{
				if (SigmaB[p][q]!=0.0)
				{
					TotalDiffDeriGinf += DifferenceDeriGinfMSF(b, ginf, eta1,eta2,rhoinf,p, q);
				}
			 }
		}
	return double((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*TotalDiffDeriGinf);
}
double TotalDerieta1(double b, double ginf,double eta1, double eta2, double rhoinf)
{

	double TotalDiffDerieta1 = 0.0;
	for( int p = 1 ; p < NUM_VALS_P ; p++)
		{
			for( int q = p+1 ; q < MAX_TENOR; q++)
			{
				if (SigmaB[p][q]!=0.0)
				{
					TotalDiffDerieta1 += DifferenceDerieta1(b, ginf, eta1,eta2,rhoinf,p, q);
				}
			 }
		}
	return double((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*TotalDiffDerieta1);
}
double TotalDerieta1MSF(double b, double ginf,double eta1, double eta2, double rhoinf)
{

	double TotalDiffDerieta1 = 0.0;
	for( int p = 1 ; p < NUM_VALS_P ; p++)
	{
		for( int q = p+1 ; q < MAX_TENOR; q++)
		{
			if (SigmaB[p][q]!=0.0)
			{
				TotalDiffDerieta1 += DifferenceDerieta1MSF(b, ginf, eta1,eta2,rhoinf,p, q);
			}
		 }
	}
	return double((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*TotalDiffDerieta1);
}
double TotalDerieta2(double b, double ginf,double eta1, double eta2, double rhoinf)
{

	double TotalDiffDerieta2 = 0.0;
	for( int p = 1 ; p < NUM_VALS_P ; p++)
		{
			for( int q = p+1 ; q < MAX_TENOR; q++)
			{
				if (SigmaB[p][q]!=0.0)
				{
					TotalDiffDerieta2 += DifferenceDerieta2(b, ginf, eta1,eta2,rhoinf,p, q);
				}
			 }
		}
	return double((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*TotalDiffDerieta2);
}
double TotalDerieta2MSF(double b, double ginf,double eta1, double eta2, double rhoinf)
{

	double TotalDiffDerieta2 = 0.0;
	for( int p = 1 ; p < NUM_VALS_P ; p++)
	{
		for( int q = p+1 ; q < MAX_TENOR; q++)
		{
			if (SigmaB[p][q]!=0.0)
			{
				TotalDiffDerieta2 += DifferenceDerieta2MSF(b, ginf, eta1,eta2,rhoinf,p, q);
			}
		 }
	}
	return double((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*TotalDiffDerieta2);
}
double TotalDerirhoinf(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	double TotalDiffDerirhoinf = 0.0;
	for( int p = 1 ; p < NUM_VALS_P ; p++)
	{
		for( int q = p+1 ; q < MAX_TENOR; q++)
		{
			if (SigmaB[p][q]!=0.0)
			{
				TotalDiffDerirhoinf += DifferenceDerirhoinf(b, ginf,eta1,eta2,rhoinf,p, q);
			}
		 }
	}
	return double((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*TotalDiffDerirhoinf);
}
double TotalDerirhoinfMSF(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	double TotalDiffDerirhoinf = 0.0;
	for( int p = 1 ; p < NUM_VALS_P ; p++)
	{
		for( int q = p+1 ; q < MAX_TENOR; q++)
		{
			if (SigmaB[p][q]!=0.0)
			{
				TotalDiffDerirhoinf += DifferenceDerirhoinfMSF(b, ginf,eta1,eta2,rhoinf,p, q);
			}
		 }
	}
	return double((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*TotalDiffDerirhoinf);
}
double RMSDeriB(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double ((1.0/CalculateRms(b, ginf,eta1,eta2,rhoinf))/2.0*TotalDeriB(b,ginf,eta1,eta2,rhoinf));
}
double RMSDeriBMSF(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double ((1.0/sqrt((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*CalculateTotalDifferenceInRmsMSF(b, ginf,eta1,eta2,rhoinf)))*TotalDeriBMSF(b,ginf,eta1,eta2,rhoinf)/2.0);
}
double RMSDeriGinf(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double ((1.0/sqrt((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*CalculateTotalDifferenceInRms(b, ginf,eta1,eta2,rhoinf)))*TotalDeriGinf(b,ginf,eta1,eta2,rhoinf)/2.0);
}
double RMSDeriGinfMSF(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double ((1.0/sqrt((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*CalculateTotalDifferenceInRmsMSF(b, ginf,eta1,eta2,rhoinf)))*TotalDeriGinfMSF(b,ginf,eta1,eta2,rhoinf)/2.0);
}
double RMSDerieta1(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double ((1.0/sqrt((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*CalculateTotalDifferenceInRms(b, ginf,eta1,eta2,rhoinf)))*TotalDerieta1(b,ginf,eta1,eta2,rhoinf)/2.0);
}
double RMSDerieta1MSF(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double ((1.0/sqrt((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*CalculateTotalDifferenceInRmsMSF(b, ginf,eta1,eta2,rhoinf)))*TotalDerieta1MSF(b,ginf,eta1,eta2,rhoinf)/2.0);
}
double RMSDerieta2(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double ((1.0/sqrt((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*CalculateTotalDifferenceInRms(b, ginf,eta1,eta2,rhoinf)))*TotalDerieta2(b,ginf,eta1,eta2,rhoinf)/2.0);
}
double RMSDerieta2MSF(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double ((1.0/sqrt((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*CalculateTotalDifferenceInRmsMSF(b, ginf,eta1,eta2,rhoinf)))*TotalDerieta2MSF(b,ginf,eta1,eta2,rhoinf)/2.0);
}
double RMSDerirhoinf(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double ((1.0/sqrt((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*CalculateTotalDifferenceInRms(b, ginf,eta1,eta2,rhoinf)))*TotalDerirhoinf(b,ginf,eta1,eta2,rhoinf)/2.0);
}
double RMSDerirhoinfMSF(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double ((1.0/sqrt((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*CalculateTotalDifferenceInRmsMSF(b, ginf,eta1,eta2,rhoinf)))*TotalDerirhoinfMSF(b,ginf,eta1,eta2,rhoinf)/2.0);
}
double MSDeriB(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double (1000.0*(2.0*CalculateRms(b,ginf,eta1,eta2,rhoinf)*RMSDeriB(b,ginf,eta1,eta2,rhoinf)*sqrt(pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4))+pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),2)*(1.0/2.0)/sqrt(pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4))*(4.0*pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),3)*RMSDeriB(b,ginf,eta1,eta2,rhoinf)+4.0*pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),3)*RMSDeriBMSF(b,ginf,eta1,eta2,rhoinf))));
}
double MSDeriB2(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double (1.0/4.0/pow((pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4)),0.75)*(4.0*pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),3)*RMSDeriB(b,ginf,eta1,eta2,rhoinf)+4.0*pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),3)*RMSDeriBMSF(b,ginf,eta1,eta2,rhoinf)));
}
double MSDeriGinf(double b, double ginf,double eta1, double eta2,  double rhoinf)
{
	return double (1000.0*(2.0*CalculateRms(b,ginf,eta1,eta2,rhoinf)*RMSDeriGinf(b,ginf,eta1,eta2,rhoinf)*sqrt(pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4))+pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),2)*(1.0/2.0)/sqrt(pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4))*(4.0*pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),3)*RMSDeriGinf(b,ginf,eta1,eta2,rhoinf)+4.0*pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),3)*RMSDeriGinfMSF(b,ginf,eta1,eta2,rhoinf))));
}
double MSDeriGinf2(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double (1.0/4.0/pow((pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4)),0.75)*(4.0*pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),3)*RMSDeriGinf(b,ginf,eta1,eta2,rhoinf)+4.0*pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),3)*RMSDeriGinfMSF(b,ginf,eta1,eta2,rhoinf)));
}
double MSDerieta1(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double (1000.0*(2.0*CalculateRms(b,ginf,eta1,eta2,rhoinf)*RMSDerieta1(b,ginf,eta1,eta2,rhoinf)*sqrt(pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4))+pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),2)*(1.0/2.0)/sqrt(pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4))*(4.0*pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),3)*RMSDerieta1(b,ginf,eta1,eta2,rhoinf)+4.0*pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),3)*RMSDerieta1MSF(b,ginf,eta1,eta2,rhoinf))));
}
double MSDerieta12(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double (1.0/4.0/pow((pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4)),0.75)*(4.0*pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),3)*RMSDerieta1(b,ginf,eta1,eta2,rhoinf)+4.0*pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),3)*RMSDerieta1MSF(b,ginf,eta1,eta2,rhoinf))); 
}
double MSDerieta2(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double (1000.0*(2.0*CalculateRms(b,ginf,eta1,eta2,rhoinf)*RMSDerieta2(b,ginf,eta1,eta2,rhoinf)*sqrt(pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4))+pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),2)*(1.0/2.0)/sqrt(pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4))*(4.0*pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),3)*RMSDerieta2(b,ginf,eta1,eta2,rhoinf)+4.0*pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),3)*RMSDerieta2MSF(b,ginf,eta1,eta2,rhoinf))));
}
double MSDerieta22(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double (1.0/4.0/pow((pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4)),0.75)*(4.0*pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),3)*RMSDerieta2(b,ginf,eta1,eta2,rhoinf)+4.0*pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),3)*RMSDerieta2MSF(b,ginf,eta1,eta2,rhoinf))); 
}
double MSDerirhoinf(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double
	(1000.0*(2.0*CalculateRms(b,ginf,eta1,eta2,rhoinf)*RMSDerirhoinf(b,ginf,eta1,eta2,rhoinf)*sqrt(pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4))+pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),2)*(1.0/2.0)/sqrt(pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4))*(4.0*pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),3)*RMSDerirhoinf(b,ginf,eta1,eta2,rhoinf)+4.0*pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),3)*RMSDerirhoinfMSF(b,ginf,eta1,eta2,rhoinf))));
}
double MSDerirhoinf2(double b, double ginf,double eta1, double eta2, double rhoinf)
{
	return double (1.0/4.0/pow((pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),4)+pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),4)),0.75)*(4.0*pow(CalculateRms(b,ginf,eta1,eta2,rhoinf),3)*RMSDerirhoinf(b,ginf,eta1,eta2,rhoinf)+4.0*pow(CalculateRmsMSF(b,ginf,eta1,eta2,rhoinf),3)*RMSDerirhoinfMSF(b,ginf,eta1,eta2,rhoinf)));  
}

double SigmaSqrDeriEta12(double eta1, double eta2, double rhoinf, int p, int q)
{
	CalculateL();
	CalculateW(p,q);
	double SigmaSqrtDeriEta1=0.0;
	for(int i=p;i<q;i++)
	{
		for(int j=p; j<q; j++)
		{
			SigmaSqrtDeriEta1 = SigmaSqrtDeriEta1 + ((W[i]*W[j]*L[i]*L[j])/(CalculateSpq( p,  q)*CalculateSpq(p,q)))*Gamma[i]*Gamma[j]*(RhoDeriEta1(eta1,eta2,rhoinf,i,j));
		}
	}
	return SigmaSqrtDeriEta1;
}
double SigmaSqrDeriEta22(double eta1, double eta2, double rhoinf, int p, int q)
{
	CalculateL();
	CalculateW(p,q);
	double SigmaSqrtDeriEta2 = 0.0;
	for(int i=p;i<q;i++)
	{
		for(int j=p; j<q; j++)
		{
			SigmaSqrtDeriEta2 = SigmaSqrtDeriEta2 + ((W[i]*W[j]*L[i]*L[j])/(CalculateSpq( p,  q)*CalculateSpq( p,  q)))*Gamma[i]*Gamma[j]*(RhoDeriEta2(eta1,eta2,rhoinf,i,j));
		}
	}
	return SigmaSqrtDeriEta2;
}
double SigmaSqrDeriRhoinf2 (double eta1, double eta2, double rhoinf, int p, int q)
{
	CalculateL();
	CalculateW(p,q);
	double SigmaSqrtDeriRhoinf = 0.0 ;
	for(int i=p;i<q;i++)
	{
		for(int j=p; j<q; j++)
		{
			SigmaSqrtDeriRhoinf = SigmaSqrtDeriRhoinf + ((W[i]*W[j]*L[i]*L[j])/(CalculateSpq( p,  q)*CalculateSpq( p,  q)))*Gamma[i]*Gamma[j]*(RhoDeriRhoinf(eta1,eta2,rhoinf, i,j));
		}
	}
	return SigmaSqrtDeriRhoinf;
}	

double SigmaDeriEta12 (double eta1,double eta2,double rhoinf,int p, int q)
{
	return double (1.0/CalculateSigma2(eta1,eta2,rhoinf, p, q)*SigmaSqrDeriEta12( eta1,eta2,rhoinf,p, q)/2.0);              
}
double SigmaDeriEta22 (double eta1,double eta2,double rhoinf,int p, int q)
{
	return double (1.0/CalculateSigma2(eta1,eta2,rhoinf, p, q)*SigmaSqrDeriEta22( eta1,eta2,rhoinf,p, q)/2.0); 
}
double SigmaDeriRhoinf2(double eta1, double eta2,double rhoinf,int p, int q)
{
	return double (1.0/CalculateSigma2( eta1,eta2,rhoinf, p,q)*SigmaSqrDeriRhoinf2 (eta1,eta2,rhoinf,p, q)/2.0);
}
double DifferenceDeriEta12(double eta1, double eta2, double rhoinf, int p, int q)
{
	return double (-2.0*(SigmaB[p][q]-CalculateSigma2(eta1,eta2,rhoinf,p,q))/(SigmaB[p][q]*SigmaB[p][q])*SigmaDeriEta12(eta1,eta2,rhoinf,p,q));
}
double DifferenceDeriEta22(double eta1, double eta2, double rhoinf, int p, int q)
{
	return double (-2.0*(SigmaB[p][q]-CalculateSigma2(eta1,eta2,rhoinf,p,q))/(SigmaB[p][q]*SigmaB[p][q])*SigmaDeriEta22(eta1,eta2,rhoinf,p,q));
}
double DifferenceDeriRhoinf2(double eta1, double eta2, double rhoinf, int p, int q)
{
	return double (-2.0*(SigmaB[p][q]-CalculateSigma2(eta1,eta2,rhoinf,p,q))/(SigmaB[p][q]*SigmaB[p][q])*SigmaDeriRhoinf2(eta1,eta2,rhoinf,p,q));
}
double TotalDeriEta12(double eta1,double eta2, double rhoinf)
{
	double TotalDiffDeriEta1 = 0.0;
	for( int p = 1 ; p < NUM_VALS_P ; p++)
		{
			for( int q = p+1 ; q < MAX_TENOR; q++)
			{
				if (SigmaB[p][q]!=0.0)
				{
					TotalDiffDeriEta1+= DifferenceDeriEta12(eta1,eta2,rhoinf, p, q);
				}
			 }
		}
	return double((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*TotalDiffDeriEta1);
}
double TotalDeriEta22(double eta1,double eta2, double rhoinf)
{
	double TotalDiffDeriEta2 = 0.0;
	for( int p = 1 ; p < NUM_VALS_P ; p++)
		{
			for( int q = p+1 ; q < MAX_TENOR; q++)
			{
				if (SigmaB[p][q]!=0.0)
				{
					TotalDiffDeriEta2+= DifferenceDeriEta22(eta1,eta2,rhoinf, p, q);
				}
			 }
		}
	return double((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*TotalDiffDeriEta2);
}
double TotalDeriRhoinf2(double eta1, double eta2, double rhoinf)
{
	double TotalDiffDeriRhoinf = 0.0;
	for( int p = 1 ; p < NUM_VALS_P ; p++)
		{
			for( int q = p+1 ; q < MAX_TENOR; q++)
			{
				if (SigmaB[p][q]!=0.0)
				{
					TotalDiffDeriRhoinf += DifferenceDeriRhoinf2(eta1,eta2,rhoinf, p, q);
				}
			 }
		}
	return double((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*TotalDiffDeriRhoinf);
}
double RMSDeriEta12(double eta1, double eta2, double rhoinf)
{
	return double ((1.0/sqrt((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*CalculateTotalDifferenceInRms2(eta1,eta2,rhoinf)))*TotalDeriEta12(eta1,eta2,rhoinf)/2.0);
}
double RMSDeriEta22(double eta1, double eta2, double rhoinf)
{
	return double ((1.0/sqrt((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*CalculateTotalDifferenceInRms2(eta1,eta2,rhoinf)))*TotalDeriEta22(eta1,eta2,rhoinf)/2.0);
}
double RMSDeriRhoinf2(double eta1, double eta2, double rhoinf)
{
	return double ((1.0/sqrt((8.0/((NUM_SWAP-2)*(NUM_SWAP-4)))*CalculateTotalDifferenceInRms2(eta1,eta2,rhoinf)))*TotalDeriRhoinf2(eta1,eta2,rhoinf)/2.0);
}


void PrintRMS(double eta1, double eta2,double rhoinf)
{
	printf ("RMS(%f, %f, %f) = %f\n", eta1, eta2,rhoinf,CalculateRms2(eta1,eta2,rhoinf));
}

class Functional : public BaseFunctional 
{
	public:
	double evaluateFG(double * x, double * g, int n);
};
double Functional::evaluateFG(double * x, double * g, int n)
{
 	double rms;
	if (l[0] == 0.0 && u[0]== 0.0)
 	{
		double eta1 = x[0];
		double eta2 = x[1];
		double rhoinf = x[2];
		rms = CalculateRms2( eta1,eta2,rhoinf);
 		g[0] = RMSDeriEta12( eta1,eta2,rhoinf);
 		g[1] = RMSDeriEta22( eta1,eta2,rhoinf);
 		g[2] = RMSDeriRhoinf2( eta1,eta2,rhoinf);
 	}
 	else
 	{
		double b = x[0];
		double ginf = x[1];
		double eta1 = x[2];
		double eta2 = x[3];
		double rhoinf = x[4];
		rms  = CalculateMS ( b, ginf,eta1,eta2,rhoinf);
		g[0] = MSDeriB ( b, ginf,eta1,eta2,rhoinf);
		g[1] = MSDeriGinf ( b, ginf,eta1,eta2,rhoinf);
		g[2] = MSDerieta1 ( b, ginf,eta1,eta2,rhoinf);
		g[3] = MSDerieta2 ( b, ginf,eta1,eta2,rhoinf);
		g[4] = MSDerirhoinf ( b, ginf,eta1,eta2,rhoinf);
 	}
	return rms;
}
double Maxerror(double b, double ginf,double eta1,double eta2,double rhoinf)
{
	double maxerror = 0.0;
	int mat = 1;
	int per = 1;
	double difference = 0.0;
	for (int p=1; p<NUM_VALS_P; p++)
	{
		for (int q=1; q< (MAX_TENOR-p);q++)
		{
			if ( SigmaB[p][q]!= 0.0)
			{
				difference = fabs((SigmaB[p][q]-CalculateSigma(b,ginf ,eta1,eta2,rhoinf,p,q)))/SigmaB[p][q];
				if (maxerror < difference)
				{
					maxerror = difference;
					mat = p/2;
					per = (q-p)/2;
				}
			}
		}
	}
	printf ("Maximum error for p=%d\t, q=%d\n", mat, per);
	return maxerror;
}

int main()
{ 
//open market data file, use default file "market_data.xml", otherwise input the file name for the market data in XML format.
	FILE *pfile; 
	string choosename;
	char filename[30],date[30],max[30];
	cout<<"________________________Joint Calibration of LIBOR Market Model:"<<std::endl;
	cout<<endl;
	loop3:
	cout<<"Data file: market_data.xml"<<endl;
 	cout<<"All Right (Ok:y, No: n)   ? ";
 	cin.clear();
 	getline(cin, choosename, '\n');
	if (choosename=="n")//open another market data file
	{
		cout<<"Data file    ? ";
		cin>>filename;
		cin.ignore();
		pfile = fopen(filename,"r");
	}
	else//open default market data file
	{
 		pfile = fopen("market_data.xml","r");	
	}
	cout<<endl;


//read the market data file	
	if (pfile == NULL)
	{
		cout << "Error: Fail to open file: "<< filename<<endl;
		goto loop3;
	}
	else
	{
		//read the generale information of market data file: date, the maximum tenor number
		ReadMark(pfile,sizeof"<date", "<date");
		fscanf(pfile, "%s", &date);
		ReadMark(pfile,sizeof "<maxnumber", "<maxnumber");
		fscanf(pfile, "%s", &max);
 		MAX_TENOR=atoi(max);
		
		//allocate memeory for vectors and arrays with size of maxnumber
		B = new double[MAX_TENOR];
		Bi = new double [MAX_TENOR];
		Gamma = new double [MAX_TENOR];
		W = new double [MAX_TENOR];
		Y = new double [MAX_TENOR];
		L = new double [MAX_TENOR];
		T = new double [MAX_TENOR];
		Shift = new double [MAX_TENOR];
		ShiftMark = new double [MAX_TENOR];
		SigmaB = new double* [MAX_TENOR];
		Smile = new double* [MAX_TENOR];
		for (int i=0; i< MAX_TENOR; i++)
		{
			SigmaB[i] = new double [MAX_TENOR+1];
			Smile[i] = new double [MAX_TENOR+1];
		}
		Mat = new int [MAX_TENOR];
		Exp = new int [MAX_TENOR];



 		//read and interpolate discount factor, or, alternatively, read swap rate and then calculate discount factor from swap rate
		ReadMark(pfile,sizeof "<dismaturity", "<dismaturity"); 
		ReadMaturity(pfile);
		fsetpos (pfile,&pos);
		if (Mat[0]!=0) //discount factor is provided
	 	{
			ReadMark(pfile,sizeof "<discountfactor", "<discountfactor");
			ReadSwapRate(pfile, Bi);
			Interpolation(MAX_TENOR,Bi);
		}
		else //discount factor is not available, read swap rate alternatively then calculate discount factor from swap rate
		{
			ReadMark(pfile,sizeof "<srmaturity", "<srmaturity");
			ReadMaturity(pfile);
			fsetpos (pfile,&pos);
			ReadMark(pfile,sizeof "<swaprate", "<swaprate");
			ReadSwapRate(pfile, B);
			Interpolation(MAX_TENOR,B);
			CalculateDisFac(MAX_TENOR);
		}
		fsetpos (pfile,&pos);

		//read and interpolate at the money(atm) caplet volatilities, if atm volatilities are not available, read out the money(otm) or in the money(itm) and from which we derive the atm volatilities.
		ReadMark(pfile,sizeof "<atmcapmaturity", "<atmcapmaturity");
		ReadMaturity(pfile);
		if (Mat[0]!=0)// atm volatilities are provided
		{
			fsetpos(pfile,&pos);
			ReadMark(pfile, sizeof"<atmcapvolatility","<atmcapvolatility");
			ReadCapVola(pfile,Gamma);
			Interpolation(MAX_TENOR-1,Gamma);
		}
		else // atm volatilities are not available, read otm or itm and then calculate atm volatilities
		{
			fsetpos(pfile,&pos);
			ReadMark(pfile,sizeof "<itmotmcapmaturity", "<itmotmcapmaturity");
			ReadMaturity(pfile);
			fsetpos(pfile,&pos);
			ReadMark(pfile, sizeof"<itmotmcapvolatility","<itmotmcapvolatility");
			ReadCapVola(pfile,Gamma);
			Interpolation(MAX_TENOR-1,Gamma);
			fsetpos(pfile,&pos);
			ReadMark(pfile,sizeof "<itmotmcapshiftmaturity", "<itmotmcapshiftmaturity");
			ReadMaturity(pfile);
			fsetpos(pfile,&pos);
			ReadMark(pfile,sizeof "<itmotmcapshift", "<itmotmcapshift");
			ReadShift(pfile);
			fsetpos(pfile,&pos);
			ReadMark(pfile,sizeof "<capshiftmaturity", "<capshiftmaturity");
			ReadMaturity(pfile);
			fsetpos(pfile,&pos);
			ReadMark(pfile,sizeof "<capshift", "<capshift");
			ReadShiftMark(pfile,"</capshift>", ShiftMark);
			fsetpos(pfile,&pos);
			ReadMark(pfile,sizeof "<capsmile", "<capsmile");
			ReadCapSmile(pfile);
			for (int i=1; i<MaxMark; i++)
			{
				Interpolation(MAX_TENOR-1,Smile[i]);
				ModifiedSmile(ShiftMark[i], Smile[i]);
			}
			CapItmOtmToAtm();// calculate atm volatilities from otm or itm volatilities 
		}
		fsetpos (pfile,&pos);

		//read and interpolate at the money(atm) swaption volatilities, if atm volatilities are not available, read out the money(otm) or in the money(itm) and from which we derive the atm volatilities.
		ReadMark(pfile,sizeof "<atmswapmaturity", "<atmswapmaturity");
		ReadMaturity(pfile);
		if (Mat[0]!=0)// atm volatilities are provided
		{
			fsetpos (pfile,&pos);
			ReadMark(pfile,sizeof "<atmswapexpiry", "<atmswapexpiry");
			ReadExpiry(pfile);
			fsetpos (pfile,&pos);
			ReadMark(pfile,sizeof "<atmswapvolatility", "<atmswapvolatility");
			ReadSwapVola(pfile);
		}
		else// atm volatilities are not available, read otm or itm and then calculate atm volatilities
		{
			fsetpos(pfile,&pos);
			ReadMark(pfile,sizeof "<itmotmswapmaturity", "<itmotmswapmaturity");
			ReadMaturity(pfile);
			fsetpos (pfile,&pos);
			ReadMark(pfile,sizeof "<itmotmswapexpiry", "<itmotmswapexpiry");
			ReadExpiry(pfile);
			fsetpos(pfile,&pos);
			ReadMark(pfile, sizeof"<itmotmswapvolatility","<itmotmswapvolatility");
			ReadSwapVola(pfile);
			fsetpos(pfile,&pos);
			ReadMark(pfile,sizeof "<itmotmswapshiftmaturity", "<itmotmswapshiftmaturity");
			ReadMaturity(pfile);
			fsetpos(pfile,&pos);
			ReadMark(pfile,sizeof "<itmotmswapshift", "<itmotmswapshift");
			ReadShift(pfile);
			fsetpos(pfile,&pos);
			ReadMark(pfile,sizeof "<swapshiftmaturity", "<swapshiftmaturity");
			ReadMaturity(pfile);
			fsetpos(pfile,&pos);
			ReadMark(pfile,sizeof "<swapshift", "<swapshift");
			ReadShiftMark(pfile,"</swapshift>", ShiftMark);
			fsetpos(pfile,&pos);
			ReadMark(pfile,sizeof "<swapsmile", "<swapsmile");
			ReadCapSmile(pfile);
			for (int i=1; i<MaxMark; i++)
			{
				Interpolation(MAX_TENOR-1,Smile[i]);
				ModifiedSmile(ShiftMark[i], Smile[i]);
			}
			SwapItmOtmToAtm();// calculate atm volatilities from otm or itm volatilities 
		}
		ReadMark(pfile,sizeof "</data", "</data");//read the ending mark of market data file and then close the file.
 		fclose(pfile);
	}

//print the general information of the input market data file and the program
	cout<<"Model: Libor Market Model (Brace-Gatarek-Musiela model) with multi-factor"<<endl;
	if (choosename=="n") cout<<"Market data:"<< filename<<endl; 
	else cout<<"Market data:	"<<"market_data.xml"<<endl;		
	cout<<"Current Date:	"<< date<<endl;
	cout<<"Tenor value (years):	0.5 "<<endl;
	cout<<"Maximum swaption maturity (in years): "<<(MAX_TENOR/2-1)<<endl;
	cout<<endl;
	cout<<"Up to maturities (by year, from 1 to "<<(MAX_TENOR/2-1)<<"):	";
	cin >> NUM_VALS_P;  // the maximum year of up to maturity of the swaption
	cin.ignore();
	NUM_VALS_P = 2*NUM_VALS_P + 1;
	DefineT();// set up the tenor structure	

//Find out the number of provided swaption volatilities within the given maximum up to maturity
	for(int i=0; i<NUM_VALS_P;i++)
	{
		for (int j=0; j<MAX_TENOR; j++)
		{
			if(SigmaB[i][j]!=0)
			{
 				NUM_SWAP++;
			}
		}
	}
 	NUM_SWAP=NUM_SWAP+1;
	

//choose parameters to be calibrated, five parameters or part of them can be calibrated, for each parameter, if it is choosen to be calibrated, the range of it should be provided, if it is not to be calibrated, a proper value for it should be provided.
	cout<<endl<<"Choose parameters to be calibrated\n"<<endl;
	cout<<"Volatility Norm g(t)=g_inf+(1-g_inf)*exp(-b*t)"<<endl;
	ChooseCaliParameter("b", 0.0, 10.00, 5.01, 0);
	ChooseCaliParameter("g_inf", 0.0, 1.00, 0.56, 1);
	cout<<"\n"<<endl;
	cout<<"Correlation structure: rho_ij(eta1,eta2,rho_inf)"<<endl;
	cout<<"\n"<<endl;
	cout<<"Choose parameters  eta1,eta2,rho_inf \n"<<endl;
	ChooseCaliParameter("eta1", 0.0, 2.00, 1.22, 2);
	ChooseCaliParameter("eta2", 0.0, 1.00, 0.001, 3);
	ChooseCaliParameter("rho_inf", 0.0, 1.00, 0.29, 4);


// Calibrate the choosen parameters,
	double minval;
	long int ll[] = {2,2,2,2,2};
	if (l[0] == 0.0 && u[0]== 0.0) // if 
	{
		Functional F;
		double x[3]={0.30,0.00,0.07};
		double l[3]={0.10,0.00,0.01};
		double u[3]={3.00,1.00,1.00};
		minval =  minimizeBFGS(x,l, u, ll, 3, F, 0.01, 0.01, 300);
		printf("RMS2(\t %f\t & %f\t & %f\t & %f\t & %f\n", x[0],x[1],x[2],CalculateRms2(x[0],x[1],x[2]), CalculateRms2(x[0],x[1],x[2]));
 		Maxerror(ParaValue[0],ParaValue[1],x[0],x[1],x[2]);

		cout<<"___________________________Calibrate eta1, eta2,rhoinf________"<<endl;
		cout<<"Date: "<<date<<endl;
		cout<<"Up to maturity: "<<(NUM_VALS_P-1)/2<<endl;
		cout<<"Number of swaption: "<<NUM_SWAP-1<<endl;
		cout<<"b: "<<ParaValue[0]<<endl;
		cout<<"ginf: "<<ParaValue[1]<<endl;
		cout<<"____________________________Result:___________________________"<<endl;
		cout<<"Eta1:     "<< x[0]<<endl;
		cout<<"Eta2:     "<< x[1]<<endl;
		cout<<"Rhoinf:   "<< x[2]<<endl;
		cout<<"Rms:      "<< CalculateRms2(x[0],x[1],x[2])<<endl;
		cout<<"Msf:      "<< CalculateRms2(x[0],x[1],x[2])<<endl;
		cout<<"Max error:"<< Maxerror(ParaValue[0],ParaValue[1],x[0],x[1],x[2])<<endl;
		}
	else 
	{
		Functional F;
		minval =  minimizeBFGS(x,l,u, ll, 5, F, 0.000001, 0.00001, 50);
		printf("\nUse the data of %s to calibrate the interest rate model\n",date);
		printf ("The results for UpToMat=%d and #swapnts= %d is \n", (NUM_VALS_P-1)/2, NUM_SWAP-1);
		printf("RMS(\t %f\t & %f\t & %f\t & %f\t & %f\t & %f\t & %f\n",x[0],x[1],x[2],x[3],x[4],CalculateRms(x[0],x[1],x[2],x[3],x[4]), CalculateRmsMSF(x[0],x[1],x[2],x[3],x[4]));
		Maxerror(x[0],x[1],x[2],x[3],x[4]);
		cout<<"___________________________Calibrate ________"<<endl;
		cout<<"Date: "<<date<<endl;
		cout<<"Up to maturity: "<<(NUM_VALS_P-1)/2<<endl;
		cout<<"Number of swaption: "<<NUM_SWAP-1<<endl;
		cout<<"____________________________Result:___________________________"<<endl;
		cout<<"b:        "<< x[0]<<endl;
		cout<<"ginf:     "<< x[1]<<endl;
		cout<<"eta1:     "<< x[2]<<endl;
		cout<<"eta2:     "<< x[3]<<endl;
		cout<<"rhoinf:   "<< x[4]<<endl;
		cout<<"RMS:      "<< CalculateRms(x[0],x[1],x[2],x[3],x[4])<<endl;
		cout<<"RMSMSF:      "<< CalculateRmsMSF(x[0],x[1],x[2],x[3],x[4])<<endl;
		cout<<"Max error:"<< Maxerror(x[0],x[1],x[2],x[3],x[4])<<endl;
	}

	//free allocations
	delete B;
	delete Bi;
	delete Gamma;
	delete W;
	delete Y;
	delete L;
	delete T;
	delete Shift;
	delete ShiftMark;
	delete Mat;
	delete Exp;
	for (int i=0; i<MAX_TENOR+1; i++)
	{
		delete SigmaB[i];
		delete Smile[i];
	}
	delete SigmaB;
	delete Smile;

	return 0;
}
