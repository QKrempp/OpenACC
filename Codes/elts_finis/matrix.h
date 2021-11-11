#ifndef MATRIX_H
#define MATRIX_H

//================================================================================
// TYPES for VECTORS and SPARSE MATRICES
//================================================================================

typedef struct Vector    Vector;
typedef struct MatrixSp  MatrixSp;

struct Vector
{
  int rows;
  double *val;
};

struct MatrixSp
{
  int rows;
  int cols;
  int nnz;
  int nnzMax;
  int *idR;
  int *idC;
  double *val;
};

//================================================================================
// FUNCTIONS for VECTORS and SPARSE MATRICES
//================================================================================

// Create

Vector*   createVector(int rows);
MatrixSp* createMatrixSp(int rows, int cols, int nnzMax);

// Destroy

void destroyVector(Vector *V);
void destroyMatrixSp(MatrixSp *M);

// Print

void printVector(Vector *V);
void printMatrixSp(MatrixSp *M);

// Sort by indices

void sortMatrixSp(MatrixSp *M);

#endif // MATRIX_H
