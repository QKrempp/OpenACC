#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

//================================================================================
// CREATE/DESTROY/PRINT VECTOR/MATRIXSP
//================================================================================

Vector* createVector(int rows){
  Vector* V = malloc(sizeof(Vector));
  V->rows = rows;
  V->val = malloc(rows * sizeof(double));
  return V;
}

MatrixSp* createMatrixSp(int rows, int cols, int nnzMax){
  MatrixSp* M = malloc(sizeof(MatrixSp));
  M->rows = rows;
  M->cols = cols;
  M->nnz = 0;
  M->nnzMax = nnzMax;
  M->idR = malloc(nnzMax * sizeof(int));
  M->idC = malloc(nnzMax * sizeof(int));
  M->val = malloc(nnzMax * sizeof(double));
  return M;
}

void destroyVector(Vector *V){
  free(V->val);
  free(V);
}

void destroyMatrixSp(MatrixSp *M){
  free(M->idR);
  free(M->idC);
  free(M->val);
  free(M);
}

void printVector(Vector *V){
  for(int i=0; i<V->rows; i++){
    printf("%f\n", V->val[i]);
  }
}

void printMatrixSp(MatrixSp *M){
  for(int n=0; n<M->nnz; n++){
    printf("(%i,%i) %f\n", M->idR[n], M->idC[n], M->val[n]);
  }
}

//================================================================================
// SORT (by indices) MATRIXSP
//================================================================================

typedef struct sortMatrixSpStr sortMatrixSpStr;
struct sortMatrixSpStr
{
  int idR;
  int idC;
  double val;
};

static int sortMatrixSpFunc(const void *a, const void *b) {
  sortMatrixSpStr *A = (sortMatrixSpStr *)a;
  sortMatrixSpStr *B = (sortMatrixSpStr *)b;
  if(A->idR < B->idR) return -2;
  if(A->idR > B->idR) return 2;
  if(A->idC < B->idC) return -1;
  if(A->idC > B->idC) return 1;
  return 0;
}

void sortMatrixSp(MatrixSp *M){
  
  sortMatrixSpStr* tmp = malloc(M->nnz * sizeof(sortMatrixSpStr));
  for(int n=0; n<M->nnz; n++){
    tmp[n].idR = M->idR[n];
    tmp[n].idC = M->idC[n];
    tmp[n].val = M->val[n];
  }
  qsort(tmp, M->nnz, sizeof(sortMatrixSpStr), sortMatrixSpFunc);
  
  for(int n=0; n<M->nnz; n++){
    M->idR[n] = tmp[n].idR;
    M->idC[n] = tmp[n].idC;
    M->val[n] = tmp[n].val;
  }
  free(tmp);
};
