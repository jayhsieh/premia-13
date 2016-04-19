#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include<malloc.h>
#include <pnl_random.h>
#include <pnl_mathtools.h>
#include <pnl_root.h>

#define TAILLE_MAX 1000

// Libors and numeraires datas
// Matrix of Markov process data

typedef struct discrete_fct
{ 
  int row, col; 
  PnlMat *grid;
  PnlMat *val;
} discrete_fct ;
// needed to use pnl_root_brent();

// params containds the prices of Digital Caplets in arreas for several strikes.
// In the case of several maturities, change PnlVect by PnlMat
typedef struct params {
  int size;
  int indexMat;		
  double cte; // this constant will represent a value of J(x*) wich will be set.
  PnlMat*V_K;
  PnlVect *K;


} params;

//Create discrete function
void Set_discrete_fct (discrete_fct *f, PnlMat *grid)
{ 

  f->row = grid->m;
  f->col = grid->n;

  f->grid = pnl_mat_new ();
  pnl_mat_clone (f->grid, grid);
  f->val = pnl_mat_create_from_double(f->row,f->col,0.0);

}

//Delete discrete function
void Delete_discrete_fct (discrete_fct *f)
{   f->col=0;
  f->row=0;
  pnl_mat_free (&(f->grid));
  pnl_mat_free (&(f->val));

}

//Set parameters
void Set_params(params *p, double cte, PnlMat *V_K, PnlVect *K)
{ 

  p->size = K->size;
  p->cte = cte;
  p->V_K=pnl_mat_copy(V_K);
  p->K=pnl_vect_copy(K);
  p->indexMat =0;
}

//Delete parameters
void Delete_params (params *p){	
  p->size=0;
  p->cte=0;
  p->indexMat =0;
  pnl_vect_free (&p->K);
  pnl_mat_free (&p->V_K);
}

void set_stateVector (PnlMat *x, PnlMat *sigma_HW, PnlVect *T){
  // each column of x is a state vector process for one state on the word w_j
  // each row of x represents a markov process 
  // sigma_HW: variance-covariance matrix of the markov process 	
  PnlMat *L;
  PnlVect *aux1,*aux2;
  int i,j,n,m;
  m = x->n;
  n = x->m;

  L=pnl_mat_copy(sigma_HW);
  pnl_mat_chol (L);
  aux1=pnl_vect_create_from_double(n,0.0);
  aux2=pnl_vect_create(n);

  for(j=0; j<m; j++){

    for(i=0; i<n; i++)
      pnl_vect_set(aux1,i,pnl_rand_normal(0));	


    pnl_mat_mult_vect_inplace(aux2,L,aux1);//aux2=x^tilde_i standard  normal correlated vector
    pnl_vect_map (aux1, T, &sqrt);// aux1_i= sqrt(T_i)= sqrt(var(x_i))
    pnl_vect_mult_vect_term (aux2, aux1);// aux2= sqrt(var(x))*x^tilde_i
    //pnl_vect_mult_double(aux2,0.01);
    pnl_mat_set_col (x,aux2,j);


  }//end j-loop
  pnl_vect_free(&aux1);
  pnl_vect_free(&aux2);
  pnl_mat_free(&L);

}

double DigitalCapletArrears (double x, void* pp){
  int indexMaturity = 2;
  // this function is the linear interpolation of the prices of the digital caplets in arreas.
  // this function give as ah result V_k(x)-J(x*) where x* is the constant set in p.
  params *p = (params*)(pp);
  int i, conter=0;
  double result;
  indexMaturity = p->indexMat;
  if(pnl_vect_get(p->K,p->size-1)<=x)
    return pnl_mat_get(p->V_K,p->size-1,indexMaturity);
  else if(pnl_vect_get(p->K,0)>x)
    return pnl_mat_get(p->V_K,0,indexMaturity);

  else{

    for(i=0; i<(p->size-1); i++)
      if(pnl_vect_get(p->K,i)<=x && pnl_vect_get(p->K,i+1)>x)
        conter=i;

    result= x-pnl_vect_get(p->K,conter);
    result= (pnl_mat_get(p->V_K,conter+1,indexMaturity)-pnl_mat_get(p->V_K,conter,indexMaturity))/(pnl_vect_get(p->K,conter+1)-pnl_vect_get(p->K,conter)) *result;
    result= result + pnl_mat_get(p->V_K,conter,indexMaturity);

    return result-p->cte;
  }
}

void sort_stateVector (PnlMat *x){
  // sorting the state matrix process will speed-up the Monte-Carlo
  pnl_mat_qsort(x,'c', 'i');
 
}

void set_grid (PnlMat *grid, PnlMat *x, int h){
  // this function set the a grid using the state matrix process
  // the step lenght in "i" is (x_imax-x_imin)/h-1;
  PnlVect *ini,*per;  
  int i,j;
  ini = pnl_vect_new ();
  per = pnl_vect_new ();

  pnl_mat_get_col (ini, x, 0);
  pnl_mat_get_col (per, x, x->n-1); 
  pnl_vect_minus_vect (per, ini);
  pnl_vect_div_double (per, h-1);

  for (i=0; i < grid->m; i++)
    for(j=0; j < h; j++)
      pnl_mat_set(grid, i, j, pnl_vect_get ( ini, i)+ j*pnl_vect_get (per,i) );

  pnl_vect_free (&ini);
  pnl_vect_free (&per);


}

void set_J (PnlVect* J_i, PnlMat *x, PnlMat *grid, PnlVect* N_i, double N0, int i){
  // J_i(x*)= N_0/m * sum(1{x_i,j> x* }/N_i(x_i,j));
  int  h, m, j, k, reference;
  double x_star;
  PnlVect *aux;

  m = x->n;
  h = grid->n;
  aux = pnl_vect_copy (N_i); 

  pnl_vect_inv_term(aux);
  pnl_vect_cumsum(aux);
  pnl_vect_mult_double(aux, N0/(double)m);
  reference = 0;

	
  for(j=0; j < h; j++){ 
    x_star= pnl_mat_get (grid, i, j);
		
    // we strongly use that x and g are sorted.
    if(x_star> pnl_mat_get (x, i, m-1))  
      pnl_vect_set (J_i, j, 0);

    else{
      for(k= reference; k < m; k++){
        if( pnl_mat_get (x, i, k) >= x_star){
          reference = k; 
          break;
        }}
      if(reference<=0)  
        pnl_vect_set (J_i, j, pnl_vect_get (aux, m-1));
      else
        pnl_vect_set (J_i, j, pnl_vect_get (aux, m-1)-pnl_vect_get (aux, reference-1)); //because x has been sorted.
    }
		
  }
  pnl_vect_free(&aux);
	


}

void Inv_HW_DigitalCapletArrears(int indexMaturity,PnlVect* L_xTilde_i, PnlVect* J_w, PnlFunc* digCapArr, double *ordre, double eps){
  /*this function serch K* such as V_i(K*)=J_i(x*)
   if J_i(x*) is bigger (smaller)  that max_V_i (min_V_i) we take K*= arg_max_V_i (K*= arg_min_V_i) */

  int i;
  double debug1,debug2;
  params* p = (params*)(digCapArr->params);
  p->indexMat = indexMaturity;
  for (i=0; i<J_w->size; i++){ 
    p->cte = pnl_vect_get(J_w,i);
    debug1=pnl_mat_get(p->V_K, 0,indexMaturity);
    debug2=pnl_mat_get(p->V_K, p->size-1,indexMaturity);
    if(p->cte >= pnl_mat_get(p->V_K, 0, indexMaturity)-eps) // beacause V_k is decreasing
      pnl_vect_set(L_xTilde_i, i, pnl_vect_get(p->K, 0));
    else if (p->cte <= pnl_mat_get(p->V_K, p->size-1,indexMaturity)+eps)
      pnl_vect_set(L_xTilde_i, i, pnl_vect_get(p->K, p->size-1));
    else
      pnl_vect_set(L_xTilde_i, i, pnl_root_brent (digCapArr,pnl_vect_get(p->K,0), pnl_vect_get(p->K,p->size-1)-eps, ordre));
    // pnl_root_brent need J_i(x*) to be ]min_V_i,max_V_i [ 
  }
	


}

void Interp_Libor (PnlVect* L_x_i , PnlVect* L_xTilde_i, PnlMat* x, PnlMat* grid, int i){
  // we are going to use natural cubic splits for interpolate the LIBOR rates.
  PnlMat *A,*A_inv;
  PnlVect *M,*M1,*y,*y1,*y2, *x_i, *grid_i;
  PnlVect *a,*b,*c,*d;
  int m,h,j,k,conter;
  double aux,aux2;
  m = x->n;
  h= grid->n;
  //initialisation of the variables
  A = pnl_mat_create_from_double(h-2,h-2, 0.0);
  A_inv = pnl_mat_new();
  y = pnl_vect_create_subvect (L_xTilde_i,  0, h-3);
  y1 = pnl_vect_create_subvect (L_xTilde_i, 1, h-2);
  y2 = pnl_vect_create_subvect (L_xTilde_i, 2, h-1);
  x_i = pnl_vect_new();
  pnl_mat_get_row (x_i, x, i);
  grid_i = pnl_vect_new();
  pnl_mat_get_row (grid_i, grid, i);
  // 
  pnl_mat_set_diag (A, 4.0, 0);
  pnl_mat_set_diag (A, 1.0, 1);
  pnl_mat_set_diag (A, 1.0, -1);
  pnl_mat_inverse_with_chol (A_inv, A);

  // y =6(h-1)*(y-2*y1+y2)
  pnl_vect_mult_double (y1, 2.0);
  pnl_vect_minus_vect (y,y1);
  pnl_vect_plus_vect (y,y2);
  pnl_vect_mult_double (y, 6.0*(h-1.0));

  // M= A^-1*y, natural conditions M_1,M_h=0;
  M = pnl_vect_create_from_double(h,0.0);
  M1= pnl_vect_new();
  pnl_mat_mult_vect_inplace (M1, A_inv, y);

  for(j=0;j<M1->size;j++)
    pnl_vect_set (M, j+1, pnl_vect_get (M1, j));



  // interpolation to all espace
  d = pnl_vect_create_subvect (L_xTilde_i, 0, h);

  b = pnl_vect_create_subvect (M, 0, h-1);
  pnl_vect_div_double (b, 2.0);


  a = pnl_vect_copy (M1);
  pnl_vect_resize (a, h-1);
  pnl_vect_minus_vect (a, pnl_vect_create_subvect (M, 0, h-1));
  pnl_vect_div_double (a, 6.0*(h-1));


  pnl_vect_resize_from_double (y, h-1, pnl_vect_get (y1, h-2));
  pnl_vect_resize_from_double (y1, h-1, pnl_vect_get (y2, h-2));
  c = pnl_vect_copy (y1);
  pnl_vect_axpby (-1.0/(h-1.0), y , 1.0/(h-1.0), c);
  pnl_vect_axpby (-(h-1.0)/6, M1, 1.0, c);
  pnl_vect_axpby (2.0*(h-1.0)/6, pnl_vect_create_subvect (M, 0, h-1), 1.0, c);


  //s(x)= s_i(x), g_i <= x < g_i+1
  //s_i(x)= a[i](x-g_i)^3 + b[i](x-g_i)^2 + c[i](x-g_i) + d[i];
  conter = 0;
  for(j=0; j<x_i->size; j++){
    aux = pnl_vect_get(x_i,j);

    if(aux >= pnl_vect_get(grid_i,h-1))
      conter = h-1;
    else{
      for(k=conter; k<h-1;k++)
        if(aux>=pnl_vect_get(grid_i,k) && aux<pnl_vect_get(grid_i,k+1))
          conter = k;

    }
    aux = aux - pnl_vect_get(grid_i,conter);
    aux2= pnl_vect_get(a,conter)*pow(aux, 3)+ pnl_vect_get(b,conter)*pow(aux, 2)+ pnl_vect_get(c,conter)*aux+  pnl_vect_get(d,conter);

    pnl_vect_set (L_x_i, j, aux2);
  }



  pnl_mat_free(&A);
  pnl_mat_free(&A_inv);
  pnl_vect_free(&M);
  pnl_vect_free(&M1);
  pnl_vect_free(&y);
  pnl_vect_free(&y1);
  pnl_vect_free(&y2);
  pnl_vect_free(&a);
  pnl_vect_free(&b);
  pnl_vect_free(&c);
  pnl_vect_free(&d);
  pnl_vect_free(&x_i);
  pnl_vect_free(&grid_i);
}

void Interp_Libor_Linear (PnlVect* L_x_i , PnlVect* L_xTilde_i, PnlMat* x, PnlMat* grid, int i){
  // linear interpolation of the LIBORs rates
  int j,k, conter=0;
  double x_ij,result;
  double debug1,debug2;

  for (j=0; j<x->n; j++){
    x_ij = pnl_mat_get(x,i,j); 
    debug1 = pnl_mat_get(grid, i, 0);
    debug2 = pnl_mat_get( grid, i, grid->n-1);
    if(debug1 >= x_ij)
      pnl_vect_set(L_x_i, j, pnl_vect_get(L_xTilde_i,0));
    else if(debug2 <= x_ij)
      pnl_vect_set(L_x_i, j, pnl_vect_get(L_xTilde_i,grid->n-1));

    else{

      for(k=0; k<(L_xTilde_i->size-1); k++)
        if(pnl_mat_get(grid, i, k)<=x_ij && pnl_mat_get(grid, i, k+1)>x_ij)
          conter=k;

      result= x_ij - pnl_mat_get(grid, i, conter);
      result= (pnl_vect_get(L_xTilde_i,conter+1)-pnl_vect_get(L_xTilde_i,conter))/(pnl_mat_get(grid, i, conter+1)-pnl_mat_get(grid, i, conter)) *result;
      result= result + pnl_vect_get(L_xTilde_i,conter);
      pnl_vect_set(L_x_i, j,result);

    }
  }

}



static void HK_iterations(PnlMat* sigma_HW, double T0, double per, double *ordre, double eps, int n, int m, int h, double D0T1, PnlFunc* digCapArr, discrete_fct* N,  discrete_fct* L_x){

  //Model's parameters
  // sigma_HW            : parameter of the HW-model representing the market ( "sigma_HW[i][j]dt= dW^iW^j", sigma_HW[i][i]=1)
  // T0                  : first HK-date
  // per                 : HK-periodicity
  // n				   : number of LIBOR rates 
	
  //Monte-Carlo parameters
  // h		       : number of points in the grid
  // m                   : number of ralisations of the state vector
  //Market Data
  // D0T1                : initial zcb prices (D0T1 = P(0,T_1)) 
  // digCapArr	       : Data of the price of the digital caplets in arreas row= expiraiton time and colums= Strikes
	
  //LIBOR functional(to be calibrated to market data) 
  // N                   : functional forms of the INVERSE of the numeraire at T[i], i.e. of 1/P(T_i,T_m), for i=0,...,m-1 
  // L_x				   : Functional LIBORS  dependant f x
	

  PnlVect *T;                          /* HK-dates */
  PnlVect *alpha;		                 /* alpha[i] = year fraction from T[i] to T[i+1] */  



  PnlMat *x;							 /* simulated state matrix process of our driven n-dimensional process, row= time and colums= states of the word*/	
  PnlMat *grid;		                 /* reference grid point at wich the functional form values will be computed*/

  PnlVect *J_i;						 /* Valeurs of J at i computed at the grid*/	
  PnlVect *N_i;						 /* functional numeraires  at i computed at the state vector process */
  PnlVect *L_xTilde_i;				 /* functional LIBORS  at i computed at the grid. */
  PnlVect *aux;						 /* auxiliar vector saving temporal information*/

  int i;





  //////////////////////////////////////////  
  // initialization of the main variables //


  T = pnl_vect_create_from_double (n, 0.0);
  for (i=0; i<n; i++)  pnl_vect_set (T, i, i * per + T0) ;
  // T[0]=T0          :  first resetting date of the LIBOR	
  // T[1],...,T[m]    :  payment dates of the LIBOR



  alpha = pnl_vect_create_from_double (n, per);




  x = pnl_mat_create(n,m);
  //simulation of a set of n-dim gaussian state vector
  set_stateVector(x, sigma_HW, T);

  //sorting the stade vector
  sort_stateVector(x);

  // reference grid points
  grid = pnl_mat_create_from_double(n,h,0.0);
  set_grid(grid, x, h);

  // setting the numeraires and the LIBORS
  Set_discrete_fct(N, x);
  Set_discrete_fct(L_x, x);

  // initialization of N[0] 
  aux = pnl_vect_create_from_double (x->n, 1.0);
  pnl_mat_set_row (N->val, aux, 0);  // N_T1 =1;
  pnl_vect_free(&aux);

 


  //////////////////////////////////////////////////
  // iterative computation of the N[0]=N_1,...,N[n-1]=N_n //
  ////////////////////////////////////////////////// 
  for(i=0; i<n; i++) 
	{
      //digCapArr.params->indexMat=i;
      /*initilisation de variables auxiliaires */

      N_i = pnl_vect_new (); 
      J_i = pnl_vect_create_from_double(h,0.0);
      L_xTilde_i = pnl_vect_create_from_double(h,0.0);
      aux = pnl_vect_create (m);



      pnl_mat_get_row (N_i, N->val, i);   //N_i is the numeraire at i, we use it to calculate the numeraire at i+1
		


      set_J(J_i, x, grid, N_i, D0T1, i); // calculating J at i computed at the grid.


      Inv_HW_DigitalCapletArrears(i,L_xTilde_i , J_i, digCapArr, ordre, eps); // calculating the LIBORS at i computed at the grid


      //Interp_Libor(aux, L_xTilde_i, x, grid ,i); // interpolating the the LIBORS, here "aux" is the Libor rate at i computed at the state vector, 
      Interp_Libor_Linear (aux , L_xTilde_i, x, grid ,i);


      pnl_mat_set_row (L_x->val, aux, i); // saving the LIBORS at i computed at the state vector.


      if(i<n-1){
		// calculating N_i+1,k= (1+ alpha_i*L_i,k)*N_i,k  
		pnl_vect_mult_double (aux, pnl_vect_get(alpha,i));
		pnl_vect_plus_double (aux, 1.0);
		pnl_vect_mult_vect_term(aux,N_i); // "aux" is now the numeraire at i+1 computed at the state vector.
		
		pnl_mat_set_row (N->val, aux, i+1); //saving the numeraire at i+1 computed at the state vector.
      }

      pnl_vect_free(&J_i); //not necesary, just in case
      pnl_vect_free(&L_xTilde_i);
      pnl_vect_free(&N_i);
      pnl_vect_free(&aux);
		

	} // end of i-loop



	



  ////////////////////////////////  
  // free the variables         //

  pnl_vect_free(&T);
  pnl_vect_free(&alpha);
  pnl_vect_free(&aux);
  pnl_mat_free (&x);
  pnl_mat_free (&grid);



  // end of: free the variables //
  ////////////////////////////////
}

double Princing_DigCapArr(discrete_fct* N,  discrete_fct* L_x ,double D0,double K, int i){
  /* princing digital caplets in arreas 
   V_i(k)= N_0 sum (1{L_x_i,j>=x*}/N_i,j)*/

  double result,aux,libor;
  int j,m;
  result = 0.0;
  m = L_x->col;
  for(j=0; j<L_x->col; j++){
    aux = D0/ pnl_mat_get(N->val,i,j);
    libor = pnl_mat_get(L_x->val,i,j);
    if( libor >= K)
      result= result+aux;
  }
  return result/m;
}



void RowFromFile(char *chaine, int numCol, PnlVect *res)
{
  int i=0;
  char delims[] = "\t";
  char *result = NULL;
  result = strtok( chaine, delims );
  while( result != NULL ) {
    pnl_vect_set(res,i, atof ( result ));

    result = strtok( NULL, delims );
    i++;
  }

}



static void CalibLiborFunctional(FILE* CovMatrix,FILE* Bond,FILE* DigCaplet, discrete_fct* N, discrete_fct* L_x,  int m, int h)
{
  int i=0,j=0;
  int LiborNumber = 3;
  int StrikeNumber = 3;
  int BondNumber = 3;
  double per=1.0/12.0;
  double T0;
  double ordre  ; 
  double eps ;
  double D0T1;

  PnlVect *RowLiborMat;
  PnlVect *RowPrices;
  PnlMat *BondPrices=NULL;

  PnlVect *aux;

  PnlMat* sigma_HW;

  PnlVect  *K, *digcap;
  PnlMat *V_K; 

  params p;
  PnlFunc digCapArr;


  char chaine[TAILLE_MAX] = "";



  ordre =0.00000001 ; 
  eps = 0.00000000000001;

  if (CovMatrix != NULL)
    {
      if(fgets(chaine, TAILLE_MAX, CovMatrix) != NULL)
        LiborNumber = (int)atof ( chaine );
      if(fgets(chaine, TAILLE_MAX, CovMatrix) != NULL)
        per = atof ( chaine );

      sigma_HW = pnl_mat_create (LiborNumber,LiborNumber);
      RowLiborMat = pnl_vect_create (LiborNumber);

      if(fgets(chaine, TAILLE_MAX, CovMatrix) != NULL)
        RowFromFile(chaine, LiborNumber, RowLiborMat);


      while (fgets(chaine, TAILLE_MAX, CovMatrix) != NULL) 
        {
          RowFromFile(chaine, LiborNumber, RowLiborMat);
          for(i=0;i<LiborNumber;i++)
            {
              pnl_mat_set(sigma_HW,j,i,pnl_vect_get(RowLiborMat,i));
            }

          j++;        
        }
      fclose(CovMatrix);
    } 


  j=0;

  if (DigCaplet != NULL)
    {
      if(fgets(chaine, TAILLE_MAX, DigCaplet) != NULL)
        StrikeNumber = (int)atof ( chaine );

      RowPrices = pnl_vect_create (LiborNumber+1);

      K = pnl_vect_create (StrikeNumber);
      V_K = pnl_mat_create (StrikeNumber,LiborNumber);

      while (fgets(chaine, TAILLE_MAX, DigCaplet) != NULL) 
        {
          RowFromFile(chaine, LiborNumber+1, RowPrices);

          pnl_vect_set(K,j,pnl_vect_get(RowPrices,0));

          for(i=1;i<=LiborNumber;i++)
            pnl_mat_set(V_K,j,i-1,pnl_vect_get(RowPrices,i));


          j++;       
        }
      fclose(DigCaplet);
    } 


  i=0;
  aux = pnl_vect_create (2);
  if (Bond != NULL)
    {
      if(fgets(chaine, TAILLE_MAX, Bond) != NULL)
        BondNumber = (int)atof ( chaine );
      BondPrices = pnl_mat_create (BondNumber,2);
      while(fgets(chaine, TAILLE_MAX, Bond) != NULL)
        {
          RowFromFile(chaine, 2, aux);
          pnl_mat_set(BondPrices,i,0,pnl_vect_get(aux,0));
          pnl_mat_set(BondPrices,i,1,pnl_vect_get(aux,1));
          i++;
        }
      fclose(Bond);
    }

  D0T1 = pnl_mat_get(BondPrices,0,1);
  T0 = pnl_mat_get(BondPrices,0,0);



  Set_params(&p, 0.0, V_K, K);


  digCapArr.function =  &DigitalCapletArrears;
  digCapArr.params = &p;



//Calibration algorihm
  HK_iterations(sigma_HW, T0,  per, &ordre, eps, LiborNumber, m, h, D0T1, &digCapArr,  N,   L_x);

  //Verification of the calibration routine repricing the
  //Digital Caplet 
  for(i=0; i<LiborNumber;i++){
    digcap = pnl_vect_create(p.size);
    for (j=0; j< p.size; j++)
      pnl_vect_set(digcap, j, Princing_DigCapArr(N, L_x,D0T1, pnl_vect_get(p.K,j), i));

    printf("\n Digital Cap market amd model prices for maturity T= %f \n",per*(i+1));
    for(j=0;j< p.size;j++)
      {
        printf("%lf\t",pnl_mat_get(V_K, j,i));
        printf("%lf\n",pnl_vect_get(digcap, j));
      }
    pnl_vect_free(&digcap);
    printf("\n");
  }


  pnl_mat_free(&sigma_HW);
  pnl_mat_free(&V_K);
  pnl_vect_free(&K);
  Delete_params(&p);


}


int main(int argc, char* argv[])
{
  //Monte-Carlo parameters
  int M;   //  : number of points in the grid
  int N_MC;   // : number of Monte-Carlo simulations
	
  discrete_fct N;  // : functional forms of the INVERSE of the numeraire at T[i], i.e. of 1/P(T_i,T_m), for i=0,...,m-1 
  discrete_fct L_x; //  : Functional LIBORS  dependant f x
	 
  N_MC = 100000;
  M = 3000;
	


  FILE* CovMatrix = NULL;
  CovMatrix = fopen("CorrelMatrix.txt", "r");

  FILE* Bond = NULL;
  Bond = fopen("Bond.txt", "r");

  FILE* DigCaplet = NULL;
  DigCaplet = fopen("DigitalCpletPrices.txt", "r");

  FILE* fichier = NULL;
  fichier = fopen("data.txt", "r");

  CalibLiborFunctional(CovMatrix,Bond,DigCaplet,  &N,  &L_x,   N_MC, M);


  return 0;
}
