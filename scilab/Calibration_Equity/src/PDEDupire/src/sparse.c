#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "sparse.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle, September 2002.
*/

int isEmpty(struct sparseMat *A){

  int i,j;

  i=0;
  while (i<A->nbCol && A->columns[i]==NULL)
    i++;
  j=0;
  while (j<A->nbLines && A->lines[j]==NULL)
    j++;
  if (i==A->nbCol && j==A->nbLines)
    return 1;
  else
    return 0;
  
}

void transpSparseMat(struct sparseMat **pt_transp_A, struct sparseMat *A){

  int i,j;
  struct cell *pt_col_A,*pt_line_transp_A;
  struct cell **pt_above;
  struct cell **sentinelles_columns,**sentinelles_lines;
  
  /* memory allocation for the matrix transp_A */
  (*pt_transp_A) = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  (*pt_transp_A)->columns = (struct cell **) malloc(A->nbLines*sizeof(struct cell *));
  sentinelles_columns = (struct cell **) malloc(A->nbLines*sizeof(struct cell *));
  (*pt_transp_A)->lines = (struct cell **) malloc(A->nbCol*sizeof(struct cell *));
  sentinelles_lines = (struct cell **) malloc(A->nbCol*sizeof(struct cell *));

  pt_above = (struct cell **) malloc(A->nbLines*sizeof(struct cell *));


  (*pt_transp_A)->nbCol = A->nbLines;
  (*pt_transp_A)->nbLines = A->nbCol; 
  for (i=0;i<A->nbLines;i++){
    sentinelles_columns[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_columns[i]->down = NULL;
    sentinelles_columns[i]->right = NULL;
  } 
  for (i=0;i<A->nbCol;i++){
    sentinelles_lines[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_lines[i]->right = NULL;
    sentinelles_lines[i]->down = NULL;
  }  

  for (i=0;i<A->nbLines;i++)
    pt_above[i] = sentinelles_columns[i];


  for (i=0;i<A->nbCol;i++){
    /* i_th line of transp_A */
    pt_line_transp_A = sentinelles_lines[i];
    pt_col_A = A->columns[i];
    while (pt_col_A != NULL){
      pt_line_transp_A->right = (struct cell *) malloc(sizeof(struct cell));
      pt_line_transp_A = pt_line_transp_A->right; 
      pt_line_transp_A->value = pt_col_A->value;
      pt_line_transp_A->line = i;
      pt_line_transp_A->col =pt_col_A->line;
      pt_line_transp_A->right = NULL;
      pt_line_transp_A->down = NULL;
      pt_above[pt_col_A->line]->down = pt_line_transp_A;
      pt_above[pt_col_A->line] = pt_above[pt_col_A->line]->down;
      pt_col_A = pt_col_A->down;
    }
    /* end of line i */
    (*pt_transp_A)->lines[i] = sentinelles_lines[i]->right;
    free(sentinelles_lines[i]);

  }
  
  for (j=0;j<A->nbLines;j++){
    (*pt_transp_A)->columns[j] = sentinelles_columns[j]->down;
    free(sentinelles_columns[j]);
  }
  
  free(sentinelles_lines);
  free(sentinelles_columns);
  free(pt_above);
      

}

void sparseTimesSparse(struct sparseMat **pt_res, struct sparseMat *A, struct sparseMat *B){
/* (*pt_res) = A*B */

  int i,j;
  double res_ij;
  struct cell *pt_lineA_i, *pt_colB_j;
  struct cell **pt_above_res;
  struct cell **sentinelles_columns,**sentinelles_lines;
  struct cell *pt_line_res;


  /* memory allocation for the matrix res */
  (*pt_res) = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  (*pt_res)->columns = (struct cell **) malloc(B->nbCol*sizeof(struct cell *));
  sentinelles_columns = (struct cell **) malloc(B->nbCol*sizeof(struct cell *));
  (*pt_res)->lines = (struct cell **) malloc(A->nbLines*sizeof(struct cell *));
  sentinelles_lines = (struct cell **) malloc(A->nbLines*sizeof(struct cell *));

  pt_above_res = (struct cell **) malloc(B->nbCol*sizeof(struct cell *));


  (*pt_res)->nbCol = B->nbCol;
  (*pt_res)->nbLines = A->nbLines; 
  for (i=0;i<B->nbCol;i++){
    sentinelles_columns[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_columns[i]->down = NULL;
    sentinelles_columns[i]->right = NULL;
  } 
  for (i=0;i<A->nbLines;i++){
    sentinelles_lines[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_lines[i]->right = NULL;
    sentinelles_lines[i]->down = NULL;
  }  

  for (i=0;i<B->nbCol;i++)
    pt_above_res[i] = sentinelles_columns[i];

  
  for (i=0;i<A->nbLines;i++){
    /* computation of the ith line of res */
    pt_line_res = sentinelles_lines[i];
    for(j=0;j<B->nbCol;j++){
      pt_lineA_i = A->lines[i];
      pt_colB_j = B->columns[j];
      res_ij = 0;
      while (pt_lineA_i != NULL && pt_colB_j != NULL){
	if (pt_lineA_i->col == pt_colB_j->line){
	  res_ij = res_ij + pt_lineA_i->value*pt_colB_j->value;
	  pt_lineA_i = pt_lineA_i->right;   
	  pt_colB_j = pt_colB_j->down;
	}
	else if (pt_lineA_i->col < pt_colB_j->line){
	  pt_lineA_i = pt_lineA_i->right;
	}
	else{
	  pt_colB_j = pt_colB_j->down;
	}
      }

      if (res_ij != 0){
	pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
	pt_line_res = pt_line_res->right;
	pt_line_res->line = i;
	pt_line_res->col = j;
	pt_line_res->value = res_ij;
	pt_line_res->down = NULL;
	pt_line_res->right = NULL;
	pt_above_res[j]->down = pt_line_res;
	pt_above_res[j] = pt_above_res[j]->down;
      }
     
      
    }
    /* end of line i */
    (*pt_res)->lines[i] = sentinelles_lines[i]->right;
    free(sentinelles_lines[i]);
    sentinelles_lines[i] = NULL;
  }
  
  for (j=0;j<B->nbCol;j++){
    (*pt_res)->columns[j] = sentinelles_columns[j]->down;
    free(sentinelles_columns[j]);
    sentinelles_columns[j] = NULL;
    pt_above_res[j] = NULL;
  }
  
  free(sentinelles_lines);
  sentinelles_lines = NULL;
  free(sentinelles_columns);
  sentinelles_columns = NULL;
  free(pt_above_res);
  pt_above_res = NULL;

}



void sparseTimesVect(double *res, struct sparseMat *A, double *B){
/* res = A*B */


  int i;
  double sum;
  struct cell *pt_lineA_i;



  for (i=0;i<A->nbLines;i++){

    pt_lineA_i = A->lines[i];
    sum = 0;
    while (pt_lineA_i != NULL){
      sum = sum + pt_lineA_i->value*B[pt_lineA_i->col];
      pt_lineA_i = pt_lineA_i->right;
    }
    res[i] = sum;

  }

}



void sparseMinusSparse(struct sparseMat **pt_res, struct sparseMat *A, struct sparseMat *B){
/* (*res) = A-B */

  int i,j;
  struct cell *pt_lineA, *pt_lineB;
  struct cell **pt_above_res;
  struct cell *pt_line_res;
  struct cell **sentinelles_columns,**sentinelles_lines;


  /* memory allocation for the matrix res */
  (*pt_res) = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  (*pt_res)->columns = (struct cell **) malloc(A->nbCol*sizeof(struct cell *));
  sentinelles_columns = (struct cell **) malloc(A->nbCol*sizeof(struct cell *));
  (*pt_res)->lines = (struct cell **) malloc(A->nbLines*sizeof(struct cell *));
  sentinelles_lines = (struct cell **) malloc(A->nbLines*sizeof(struct cell *));

  (*pt_res)->nbCol = A->nbCol;
  (*pt_res)->nbLines = A->nbLines;
  
  /* memory allocation and initialization of pt_above_res */
  pt_above_res = (struct cell **) malloc(A->nbCol*sizeof(struct cell *));


  for (i=0;i<A->nbCol;i++){
    sentinelles_columns[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_columns[i]->down = NULL;
    sentinelles_columns[i]->right = NULL;    
    pt_above_res[i] = sentinelles_columns[i];
  } 
  for (i=0;i<A->nbLines;i++){
    sentinelles_lines[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_lines[i]->right = NULL;
    sentinelles_lines[i]->down = NULL;
  }  
  
  
  
  for (i=0;i<A->nbLines;i++){
    
    /* computation of the ith line of res */
    pt_line_res = sentinelles_lines[i];
    
    
    pt_lineA = A->lines[i];
    pt_lineB = B->lines[i];
    while (pt_lineA != NULL && pt_lineB != NULL){
      
      if (pt_lineA->col < pt_lineB->col){

	/* A[i][pt_lineA->col] != 0 and B[i][pt_lineA->col] == 0 */
	/* so res[i][pt_lineA->col] =  A[i][pt_lineA->col] */
	pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
	pt_line_res = pt_line_res->right;
	pt_line_res->value = pt_lineA->value;
	pt_line_res->line = i;
	pt_line_res->col = pt_lineA->col;
	pt_line_res->down = NULL;
	pt_line_res->right = NULL;
	pt_above_res[pt_lineA->col]->down = pt_line_res;
	pt_above_res[pt_lineA->col] = pt_above_res[pt_lineA->col]->down;
	pt_lineA = pt_lineA->right;

      }

      else if (pt_lineA->col > pt_lineB->col){
	/* A[i][pt_lineB->col] == 0 and B[i][pt_lineB->col] != 0 */
	/* so res[i][pt_lineB->col] =  -B[i][pt_lineB->col] */

	pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
	pt_line_res = pt_line_res->right;
	pt_line_res->value = -pt_lineB->value;
	pt_line_res->line = i;
	pt_line_res->col = pt_lineB->col;
	pt_line_res->down = NULL;
	pt_line_res->right = NULL;
	pt_above_res[pt_lineB->col]->down = pt_line_res;
	pt_above_res[pt_lineB->col] = pt_above_res[pt_lineB->col]->down;
	pt_lineB = pt_lineB->right;

      }

      else{ /*if (pt_lineA->col == pt_lineB->col) */
	
	/* A[i][pt_lineB->col] != 0 and B[i][pt_lineB->col] != 0 */
	/* so res[i][pt_lineA->col] =  A[i][pt_lineA->col]-B[i][pt_lineA->col] */
	if (pt_lineA->value-pt_lineB->value != 0){
	  pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
	  pt_line_res = pt_line_res->right;
	  pt_line_res->value = pt_lineA->value-pt_lineB->value;
	  pt_line_res->line = i;
	  pt_line_res->col = pt_lineA->col;
	  pt_line_res->down = NULL;
	  pt_line_res->right = NULL;
	  pt_above_res[pt_lineA->col]->down = pt_line_res;
	  pt_above_res[pt_lineA->col] = pt_above_res[pt_lineA->col]->down;
	}
	pt_lineA = pt_lineA->right;
	pt_lineB = pt_lineB->right;
      }

    }
    if (pt_lineA != NULL && pt_lineB == NULL){
      /* we end up the line with the elements of A */
      while (pt_lineA != NULL){
	pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
	pt_line_res = pt_line_res->right;
	pt_line_res->value = pt_lineA->value;
	pt_line_res->line = i;
	pt_line_res->col = pt_lineA->col;
	pt_line_res->down = NULL;
	pt_line_res->right = NULL;
	pt_above_res[pt_lineA->col]->down = pt_line_res;
	pt_above_res[pt_lineA->col] = pt_above_res[pt_lineA->col]->down;
	pt_lineA = pt_lineA->right;
      }
      
    }
    else if (pt_lineA == NULL && pt_lineB != NULL){
      /* we end up the line with the elements of B */
      while (pt_lineB != NULL){
	pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
	pt_line_res = pt_line_res->right;
	pt_line_res->value = -pt_lineB->value;
	pt_line_res->line = i;
	pt_line_res->col = pt_lineB->col;
	pt_line_res->down = NULL;
	pt_line_res->right = NULL;
	pt_above_res[pt_lineB->col]->down = pt_line_res;
	pt_above_res[pt_lineB->col] = pt_above_res[pt_lineB->col]->down;
	pt_lineB = pt_lineB->right;
      }
      
    }
    
    /* end of line i */
    (*pt_res)->lines[i] = sentinelles_lines[i]->right;
    free(sentinelles_lines[i]);
    sentinelles_lines[i] = NULL;

  }
  for (i=0;i<A->nbCol;i++){
    (*pt_res)->columns[i] = sentinelles_columns[i]->down;
    free(sentinelles_columns[i]);
    sentinelles_columns[i] = NULL;
    pt_above_res[i] = NULL;
  }
  
  free(sentinelles_columns);
  sentinelles_columns = NULL;
  free(sentinelles_lines);
  sentinelles_lines = NULL;
  free(pt_above_res);
  pt_above_res = NULL;
}


void sparsePlusSparse(struct sparseMat **pt_res, struct sparseMat *A, struct sparseMat *B){
/* (*res) = A+B */

  int i,j;
  struct cell *pt_lineA, *pt_lineB;
  struct cell **pt_above_res;
  struct cell *pt_line_res;
  struct cell **sentinelles_columns,**sentinelles_lines;

  /* memory allocation for the matrix res */
  (*pt_res) = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  (*pt_res)->columns = (struct cell **) malloc(A->nbCol*sizeof(struct cell *));
  sentinelles_columns = (struct cell **) malloc(A->nbCol*sizeof(struct cell *));
  (*pt_res)->lines = (struct cell **) malloc(A->nbLines*sizeof(struct cell *));
  sentinelles_lines = (struct cell **) malloc(A->nbLines*sizeof(struct cell *));

  (*pt_res)->nbCol = A->nbCol;
  (*pt_res)->nbLines = A->nbLines;
  
  /* memory allocation and initialization of pt_above_res */
  pt_above_res = (struct cell **) malloc(A->nbCol*sizeof(struct cell *));


  for (i=0;i<A->nbCol;i++){
    sentinelles_columns[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_columns[i]->down = NULL;
    sentinelles_columns[i]->right = NULL;
    pt_above_res[i] = sentinelles_columns[i];
  } 
  for (i=0;i<A->nbLines;i++){
    sentinelles_lines[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_lines[i]->right = NULL;
    sentinelles_lines[i]->down = NULL;
  }  
  
  
  
  for (i=0;i<A->nbLines;i++){
    
    /* computation of the ith line of res */
    pt_line_res = sentinelles_lines[i];
    
    
    pt_lineA = A->lines[i];
    pt_lineB = B->lines[i];
    while (pt_lineA != NULL && pt_lineB != NULL){
      
      if (pt_lineA->col < pt_lineB->col){

	/* A[i][pt_lineA->col] != 0 and B[i][pt_lineA->col] == 0 */
	/* so res[i][pt_lineA->col] =  A[i][pt_lineA->col] */

	pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
	pt_line_res = pt_line_res->right;
	pt_line_res->value = pt_lineA->value;
	pt_line_res->line = i;
	pt_line_res->col = pt_lineA->col;
	pt_line_res->down = NULL;
	pt_line_res->right = NULL;
	pt_above_res[pt_lineA->col]->down = pt_line_res;
	pt_above_res[pt_lineA->col] = pt_above_res[pt_lineA->col]->down;
	pt_lineA = pt_lineA->right;

      }

      else if (pt_lineA->col > pt_lineB->col){

	/* A[i][pt_lineB->col] == 0 and B[i][pt_lineB->col] != 0 */
	/* so res[i][pt_lineB->col] =  B[i][pt_lineB->col] */

	pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
	pt_line_res = pt_line_res->right;
	pt_line_res->value = pt_lineB->value;
	pt_line_res->line = i;
	pt_line_res->col = pt_lineB->col;
	pt_line_res->down = NULL;
	pt_line_res->right = NULL;
	pt_above_res[pt_lineB->col]->down = pt_line_res;
	pt_above_res[pt_lineB->col] = pt_above_res[pt_lineB->col]->down;
	pt_lineB = pt_lineB->right;

      }

      else{ /*if (pt_lineA->col == pt_lineB->col) */
	
	/* A[i][pt_lineB->col] != 0 and B[i][pt_lineB->col] != 0 */
	/* so res[i][pt_lineA->col] =  A[i][pt_lineA->col]+B[i][pt_lineA->col] */

	pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
	pt_line_res = pt_line_res->right;
	pt_line_res->value = pt_lineB->value+pt_lineA->value;
	pt_line_res->line = i;
	pt_line_res->col = pt_lineA->col;
	pt_line_res->down = NULL;
	pt_line_res->right = NULL;
	pt_above_res[pt_lineA->col]->down = pt_line_res;
	pt_above_res[pt_lineA->col] = pt_above_res[pt_lineA->col]->down;
	pt_lineA = pt_lineA->right;
	pt_lineB = pt_lineB->right;

      }

    }

    if (pt_lineA != NULL && pt_lineB == NULL){
      
      /* we end up the line with the elements of A */
      while (pt_lineA != NULL){
	pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
	pt_line_res = pt_line_res->right;
	pt_line_res->value = pt_lineA->value;
	pt_line_res->line = i;
	pt_line_res->col = pt_lineA->col;
	pt_line_res->down = NULL;
	pt_line_res->right = NULL;
	pt_above_res[pt_lineA->col]->down = pt_line_res;
	pt_above_res[pt_lineA->col] = pt_above_res[pt_lineA->col]->down;
	pt_lineA = pt_lineA->right;
      }
      
    }
    else if (pt_lineA == NULL && pt_lineB != NULL){
      
      /* we end up the line with the elements of B */
      while (pt_lineB != NULL){
	pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
	pt_line_res = pt_line_res->right;
	pt_line_res->value = pt_lineB->value;
	pt_line_res->line = i;
	pt_line_res->col = pt_lineB->col;
	pt_line_res->down = NULL;
	pt_line_res->right = NULL;
	pt_above_res[pt_lineB->col]->down = pt_line_res;
	pt_above_res[pt_lineB->col] = pt_above_res[pt_lineB->col]->down;
	pt_lineB = pt_lineB->right;
      }
      
    }
    
    
    /* end of line i */
    (*pt_res)->lines[i] = sentinelles_lines[i]->right;;
    free(sentinelles_lines[i]);
    sentinelles_lines[i] = NULL;
    
  }
  
  for (i=0;i<A->nbCol;i++){
    (*pt_res)->columns[i] = sentinelles_columns[i]->down;
    free(sentinelles_columns[i]);
    sentinelles_columns[i] = NULL;
    pt_above_res[i] = NULL;
  }
  
  free(sentinelles_lines);
  sentinelles_lines = NULL;
  free(sentinelles_columns);
  sentinelles_columns = NULL;
  free(pt_above_res);
  pt_above_res = NULL;
  
}






void affectSparse(struct sparseMat **pt_res, struct sparseMat *A){
  /* (*pt_res) = A */


  int i;
  struct cell *pt_lineA;
  struct cell **pt_above_res;
  struct cell *pt_line_res;
  struct cell **sentinelles_columns,**sentinelles_lines;

  /* memory allocation for the sentinelles */
  sentinelles_columns = (struct cell **) malloc(A->nbCol*sizeof(struct cell *));
  sentinelles_lines = (struct cell **) malloc(A->nbLines*sizeof(struct cell *));
  
  /* memory allocation and initialization of pt_above_res */
  pt_above_res = (struct cell **) malloc(A->nbCol*sizeof(struct cell *));

  /* init of the sentinelles */
  for (i=0;i<A->nbCol;i++){
    sentinelles_columns[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_columns[i]->down = NULL;
    sentinelles_columns[i]->right = NULL;
    pt_above_res[i] = sentinelles_columns[i];
  } 
  for (i=0;i<A->nbLines;i++){
    sentinelles_lines[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_lines[i]->right = NULL;
    sentinelles_lines[i]->down = NULL;
  }  
  
  
  
  for (i=0;i<A->nbLines;i++){
    
    /* computation of the ith line of res */
    pt_line_res = sentinelles_lines[i];
    
    
    pt_lineA = A->lines[i];
    while (pt_lineA != NULL){
      pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
      pt_line_res = pt_line_res->right;
      pt_line_res->value = pt_lineA->value;
      pt_line_res->line = i;
      pt_line_res->col = pt_lineA->col;
      pt_line_res->down = NULL;
      pt_line_res->right = NULL;
      pt_above_res[pt_lineA->col]->down = pt_line_res;
      pt_above_res[pt_lineA->col] = pt_above_res[pt_lineA->col]->down;
      pt_lineA = pt_lineA->right;
    }

    /* end of line i */
    (*pt_res)->lines[i] = sentinelles_lines[i]->right;
    free(sentinelles_lines[i]);
    sentinelles_lines[i] = NULL;

  }

  for (i=0;i<A->nbCol;i++){
    (*pt_res)->columns[i] = sentinelles_columns[i]->down;
    free(sentinelles_columns[i]);
    sentinelles_columns[i] = NULL;
    pt_above_res[i] = NULL;
  }
  
  free(sentinelles_lines);
  sentinelles_lines = NULL;
  free(sentinelles_columns);
  sentinelles_columns = NULL;
  free(pt_above_res);
  pt_above_res = NULL;

}


void regPlusSparse(double **sum, struct sparseMat *A, double **B){
/* computes sum = A+B */

  int col,line,k;
  struct cell *pt_lineA;
  
  for (line=0;line<A->nbLines;line++){
    
    pt_lineA = A->lines[line];
    col = 0;
    while (pt_lineA!=NULL){
      for (k=col;k<pt_lineA->col;k++)
	sum[line][k] = B[line][k];
      sum[line][pt_lineA->col] = pt_lineA->value + B[line][pt_lineA->col];
      col = pt_lineA->col + 1;
      pt_lineA = pt_lineA->right;
    }
    for (k=col;k<A->nbCol;k++)
      sum[line][k] = B[line][k];

  }
    
}
  

void regTimesSparse(struct sparseMat **pt_res, double **A, struct sparseMat *B, int nbLinesA){
/* computes (*pt_res) = A*B */
 
  int i,j;
  double res_ij;
  struct cell *pt_colB_j;
  struct cell **pt_above_res; 
  struct cell **sentinelles_columns,**sentinelles_lines;
  struct cell *pt_line_res;


  /* memory allocation for the matrix res */
  (*pt_res) = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  (*pt_res)->columns = (struct cell **) malloc(B->nbCol*sizeof(struct cell *));
  sentinelles_columns = (struct cell **) malloc(B->nbCol*sizeof(struct cell *));
  (*pt_res)->lines = (struct cell **) malloc(nbLinesA*sizeof(struct cell *));
  sentinelles_lines = (struct cell **) malloc(nbLinesA*sizeof(struct cell *));
  (*pt_res)->nbCol = B->nbCol;
  (*pt_res)->nbLines = nbLinesA;
  for (i=0;i<B->nbCol;i++){
    sentinelles_columns[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_columns[i]->down = NULL;    
    sentinelles_columns[i]->right = NULL;    
  }
  for (i=0;i<nbLinesA;i++){
    sentinelles_lines[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_lines[i]->right = NULL;  
    sentinelles_lines[i]->down = NULL;  

  }

  
  /* memory allocation and initialization of pt_above_res */
  pt_above_res = (struct cell **) malloc(B->nbCol*sizeof(struct cell *));

  for (i=0;i<B->nbCol;i++)
    pt_above_res[i] = sentinelles_columns[i];


  for (i=0;i<nbLinesA;i++){
    /* computation of the ith line of res */
    pt_line_res = sentinelles_lines[i];
    for(j=0;j<B->nbCol;j++){
      pt_colB_j = B->columns[j];
      res_ij = 0;
      while (pt_colB_j != NULL){
	res_ij = res_ij + A[i][pt_colB_j->line]*pt_colB_j->value;
	pt_colB_j = pt_colB_j->down;
      }
      if (res_ij != 0){
	pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
	pt_line_res = pt_line_res->right;
	pt_line_res->line = i;
	pt_line_res->col = j;
	pt_line_res->value = res_ij;
	pt_line_res->down = NULL;
	pt_line_res->right = NULL;
	pt_above_res[j]->down = pt_line_res;
	pt_above_res[j] = pt_above_res[j]->down;
      }

    }
    /* end of line i */
    (*pt_res)->lines[i] = sentinelles_lines[i]->right;
    free(sentinelles_lines[i]);
    sentinelles_lines[i] = NULL;
  }
 
  for (j=0;j<B->nbCol;j++){
    (*pt_res)->columns[j] = sentinelles_columns[j]->down;
    free(sentinelles_columns[j]);
    sentinelles_columns[j] = NULL;
    pt_above_res[j] = NULL;
  }
  
  free(sentinelles_lines);
  sentinelles_lines = NULL;
  free(sentinelles_columns);
  sentinelles_columns = NULL;
  free(pt_above_res);
  pt_above_res = NULL;

}



void sparseMatMult(struct sparseMat **pt_res, struct sparseMat *mat){
/* (*res) = transpose(mat)*mat */
 
  int i,j;
  double res_ij;
  struct cell *pt_col_i, *pt_col_j;
  struct cell **pt_above_res;
  struct cell **sentinelles_columns,**sentinelles_lines;
  struct cell *pt_line_res;


  (*pt_res) = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  (*pt_res)->columns = (struct cell **) malloc(mat->nbCol*sizeof(struct cell *));
  sentinelles_columns = (struct cell **) malloc(mat->nbCol*sizeof(struct cell *));
  (*pt_res)->lines = (struct cell **) malloc(mat->nbCol*sizeof(struct cell *));
  sentinelles_lines = (struct cell **) malloc(mat->nbCol*sizeof(struct cell *));  
  (*pt_res)->nbCol = mat->nbCol;
  (*pt_res)->nbLines = mat->nbCol;
  for (i=0;i<mat->nbCol;i++){
    sentinelles_columns[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_columns[i]->down = NULL;
    sentinelles_columns[i]->right = NULL;    
  }

  for (i=0;i<mat->nbCol;i++){
    sentinelles_lines[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_lines[i]->right = NULL;
    sentinelles_lines[i]->down = NULL;  
  }

  pt_above_res = (struct cell **) malloc(mat->nbCol*sizeof(struct cell *));

  for (i=0;i<mat->nbCol;i++)
    pt_above_res[i] = sentinelles_columns[i];

  for (i=0;i<mat->nbCol;i++){
    /* computation of the ith line of res */
    pt_line_res = sentinelles_lines[i];
 
    for(j=0;j<mat->nbCol;j++){
      
      pt_col_i = mat->columns[i];
      pt_col_j = mat->columns[j];
      res_ij = 0;
      while (pt_col_i != NULL && pt_col_j != NULL){
	if (pt_col_i->line == pt_col_j->line){
	  res_ij = res_ij + pt_col_i->value*pt_col_j->value;
	  pt_col_i = pt_col_i->down;
	  pt_col_j = pt_col_j->down;
	}
	else if (pt_col_i->line < pt_col_j->line)
	  pt_col_i = pt_col_i->down;
	else
	  pt_col_j = pt_col_j->down;
      }
      
      if (res_ij != 0){
	pt_line_res->right = (struct cell *) malloc(sizeof(struct cell));
	pt_line_res = pt_line_res->right;
	pt_line_res->line = i;
	pt_line_res->col = j;
	pt_line_res->value = res_ij;
	pt_line_res->down = NULL;
	pt_line_res->right = NULL;
	pt_above_res[j]->down = pt_line_res;
	pt_above_res[j] = pt_above_res[j]->down;
      }
     
      
    }
    /* end of line i */
    (*pt_res)->lines[i] = sentinelles_lines[i]->right;
    free(sentinelles_lines[i]);
    sentinelles_lines[i] = NULL;
  }
  
  for (j=0;j<mat->nbCol;j++){
    (*pt_res)->columns[j] = sentinelles_columns[j]->down;
    free(sentinelles_columns[j]);
    sentinelles_columns[j] = NULL;
    pt_above_res[j] = NULL;
  }
  
  free(sentinelles_lines);
  sentinelles_lines = NULL;
  free(sentinelles_columns);
  sentinelles_columns = NULL;
  free(pt_above_res);
  pt_above_res = NULL;

} 



void desallocSparseMat_cells(struct sparseMat *mat){
/* frees the cells of the sparse matrix "mat" */

  struct cell *pt_prev, *pt;
  int line,col;

  for (line=mat->nbLines-1;line>=0;line--){

    if (mat->lines[line] != NULL){
      pt_prev = mat->lines[line];
      pt = pt_prev;
      while (pt != NULL){
	pt = pt->right;
	free(pt_prev);
	pt_prev = pt;
      }
      mat->lines[line] = NULL;
    }

  }
    
  for (col=0;col<mat->nbCol;col++)
    mat->columns[col] = NULL;

}



void printSparseMat(struct sparseMat *mat){
/* prints the sparse matrix "mat" */

  int i,j,k;
  struct cell *pt_cell;

  for (i=0;i<mat->nbLines;i++){
    //printf("line=%d\n",i);
    k=0;
    pt_cell = mat->lines[i];
    /* if line i is empty */
    if (pt_cell == NULL){
      for (j=0;j<mat->nbCol;j++)
	printf("%lf  ",(double)0);
    }
    else{
      while (pt_cell != NULL){
	for (j=k;j<pt_cell->col;j++)
	  printf("%lf ",(double)0);
	printf("%lf  ",pt_cell->value);
	if (pt_cell->right == NULL)
	  for (j=pt_cell->col+1;j<mat->nbCol;j++)
	    printf("%lf ",(double)0);	  
	k = pt_cell->col+1;     
	pt_cell = pt_cell->right;
      }
    }
    printf("\n");
  }
  
}


void testSparseMat(struct sparseMat **pt_mat){
/*builds a test sparse matrix (*pt_mat) */

  int i;
  struct cell *pt_cell;

  (*pt_mat) = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  (*pt_mat)->columns = (struct cell **) malloc(10*sizeof(struct cell *));
  (*pt_mat)->lines = (struct cell **) malloc(2*sizeof(struct cell *));
  for (i=0;i<10;i++)
    (*pt_mat)->columns[i] = (struct cell *) malloc(sizeof(struct cell));
  for (i=0;i<2;i++)
    (*pt_mat)->lines[i] = (struct cell *) malloc(sizeof(struct cell));
  (*pt_mat)->nbLines = 2;
  (*pt_mat)->nbCol = 10;
 
  (*pt_mat)->columns[0] = NULL;
  pt_cell = (*pt_mat)->lines[0];
  pt_cell->value = 1;
  pt_cell->col = 1;
  pt_cell->line = 0;
  pt_cell->down = NULL;
  (*pt_mat)->columns[1] = pt_cell;
  (*pt_mat)->columns[2] = NULL;
  (*pt_mat)->columns[4] = NULL;  
  pt_cell->right = (struct cell *) malloc(sizeof(struct cell));
  pt_cell = pt_cell->right;
  pt_cell->value = 2;
  pt_cell->col = 5;
  pt_cell->line = 0;
  pt_cell->down = NULL;
  (*pt_mat)->columns[5] = pt_cell;
  pt_cell->right = (struct cell *) malloc(sizeof(struct cell));
  pt_cell = pt_cell->right;
  pt_cell->value = 3;
  pt_cell->col = 6;
  pt_cell->line = 0;
  (*pt_mat)->columns[6] = pt_cell;
  (*pt_mat)->columns[7] = NULL;  
  pt_cell->right = (struct cell *) malloc(sizeof(struct cell));
  pt_cell = pt_cell->right;
  pt_cell->value = 4;
  pt_cell->col = 9;
  pt_cell->line = 0;
  (*pt_mat)->columns[9] = pt_cell;
  pt_cell->down = NULL;
  pt_cell->right = NULL;

  pt_cell = (*pt_mat)->lines[1];
  pt_cell->value = 5;
  pt_cell->col = 3;
  pt_cell->line = 1;
  pt_cell->down = NULL;
  (*pt_mat)->columns[3] = pt_cell;
  pt_cell->right = (struct cell *) malloc(sizeof(struct cell));
  pt_cell = pt_cell->right;
  pt_cell->value = 6;
  pt_cell->col = 6;
  pt_cell->line = 1;
  pt_cell->down = NULL;
  (*pt_mat)->columns[6]->down = pt_cell;
  pt_cell->right = (struct cell *) malloc(sizeof(struct cell));
  pt_cell = pt_cell->right;
  pt_cell->value = 7;
  pt_cell->col = 8;
  pt_cell->line = 1;
  pt_cell->down = NULL;
  (*pt_mat)->columns[8] = pt_cell;
  pt_cell->right = NULL;
  

}

void testSparseMat2(struct sparseMat **pt_mat){

  int i;
  struct cell *pt_cell;

  (*pt_mat) = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  (*pt_mat)->columns = (struct cell **) malloc(2*sizeof(struct cell *));
  (*pt_mat)->lines = (struct cell **) malloc(2*sizeof(struct cell *));
  for (i=0;i<2;i++)
    (*pt_mat)->columns[i] = (struct cell *) malloc(sizeof(struct cell));
  for (i=0;i<2;i++)
    (*pt_mat)->lines[i] = (struct cell *) malloc(sizeof(struct cell));
  (*pt_mat)->nbLines = 2;
  (*pt_mat)->nbCol = 2;
 
  pt_cell = (*pt_mat)->lines[0];
  pt_cell->value = 1;
  pt_cell->col = 1;
  pt_cell->line = 0;
  pt_cell->down = NULL;
  (*pt_mat)->columns[1] = pt_cell;
  pt_cell->right = NULL;

  pt_cell = (*pt_mat)->lines[1];
  pt_cell->value = 2;
  pt_cell->col = 0;
  pt_cell->line = 1;
  pt_cell->down = NULL;
  (*pt_mat)->columns[0] = pt_cell;
  pt_cell->right = NULL;

}


void tridiagToSparse(struct tridiag *T, struct sparseMat **pt_res){
/* (*pt_res) = T */

  /*declarations*/
  int col,line;
  struct cell **above,**sentinelles_columns;
  struct cell *newCell;

  /*memory allocation*/
  above = (struct cell **)malloc(T->size*sizeof(struct cell *));
  sentinelles_columns = (struct cell **)malloc(T->size*sizeof(struct cell *));
 
  (*pt_res) = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  (*pt_res)->columns = (struct cell **) malloc(T->size*sizeof(struct cell *));
  (*pt_res)->lines = (struct cell **) malloc((T->size-2)*sizeof(struct cell *));
  (*pt_res)->nbCol = T->size;
  (*pt_res)->nbLines = T->size-2;

  for (col=0;col<T->size;col++){
    sentinelles_columns[col] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_columns[col]->down = NULL;
    sentinelles_columns[col]->right = NULL;    
  }

  for (col=0;col<T->size;col++)
     above[col] = sentinelles_columns[col];
 
  for (line=0;line<T->size-2;line++){
    (*pt_res)->lines[line] = (struct cell *)malloc(sizeof(struct cell));
    newCell = (*pt_res)->lines[line];
    newCell->value = T->subdiag[line];
    newCell->line = line;
    newCell->col = line;
    newCell->down = NULL;
    above[line]->down = newCell;
    above[line] = above[line]->down;
    
    newCell->right = (struct cell *)malloc(sizeof(struct cell));
    newCell = newCell->right;
    newCell->value = T->diag[line+1];
    newCell->line = line;
    newCell->col = line+1;
    newCell->down = NULL;    
    above[line+1]->down = newCell;
    above[line+1] = above[line+1]->down;
    
    newCell->right = (struct cell *)malloc(sizeof(struct cell));
    newCell = newCell->right;
    newCell->value = T->updiag[line+1];
    newCell->line = line;
    newCell->col = line+2;
    newCell->right = NULL;
    above[line+2]->down = newCell;
    above[line+2] = above[line+2]->down;
 
  }

  for (col=0;col<T->size;col++){
    (*pt_res)->columns[col] = sentinelles_columns[col];
    free(sentinelles_columns[col]);
    sentinelles_columns[col] = NULL;
    above[col] = NULL;
  }

  free(sentinelles_columns);
  sentinelles_columns = NULL;
  free(above);
  above = NULL;

}


void scalarTimesSparse(double a, struct sparseMat **M){
  /* computes (*M) = a*(*M)   */

  int i;
  struct cell *pt_cell;

  for (i=0;i<(*M)->nbLines;i++){
    /* for each line of M */
    pt_cell = (*M)->lines[i];
    while (pt_cell != NULL){
      pt_cell->value = a*pt_cell->value;
      pt_cell = pt_cell->right;
    }
  }

}
