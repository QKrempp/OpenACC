// icc -std=c99 matrix.c mesh.c main.c -o myfem

#include <math.h>
#include <stdio.h>
#include <accel.h>

#include "matrix.h"
#include "mesh.h"


int main(int argc, char* argv[])
{
  
  // ------------------------------------------------------------------------------
  // 1. Build MESH
  // ------------------------------------------------------------------------------
  
  printf("== read mesh\n");
  
  Mesh *mesh = readMsh("square.msh");
  
  // ------------------------------------------------------------------------------
  // 2. Build FIELDS
  // ------------------------------------------------------------------------------
  
  printf("== build fields\n");
  
  const int numDof = mesh->vertices->num;
  Vector *solNum = createVector(numDof);
  Vector *solRef = createVector(numDof);
  Vector *f      = createVector(numDof);
  
  for(int i=0; i<numDof; ++i){
    double x = mesh->vertices->x[i];
    double y = mesh->vertices->y[i];
    solNum->val[i] = 0.;
    solRef->val[i] = cos(M_PI*x)*cos(2*M_PI*y);
    f->val[i]      = (1+5*M_PI*M_PI)*solRef->val[i];
  }
  
  saveToMsh(solNum->val, mesh, "solInit", "solInit.msh");
  saveToMsh(solRef->val, mesh, "solRef",  "solRef.msh");
  
  // ------------------------------------------------------------------------------
  // 3. Build LOCAL & GLOBAL Mass/Stiffness MATRICES
  // ------------------------------------------------------------------------------
  
  printf("== build mass/stiffness matrices\n");
  
  MatrixSp *M = createMatrixSp(numDof, numDof, mesh->elements->num * 9);
  MatrixSp *K = createMatrixSp(numDof, numDof, mesh->elements->num * 9);
  
  printf("   ... build local matrices\n");
  
  for(int e=0; e<mesh->elements->num; ++e){
    
    // Get infos on nodes
    int s[3];
    double x[3], y[3];
    for(int n=0; n<3; ++n){
      s[n] = mesh->elements->vertices[e][n];
      x[n] = mesh->vertices->x[s[n]];
      y[n] = mesh->vertices->y[s[n]];
    }
    
    // Compute local matrices
    double tmp[3][2];
    tmp[0][0] = y[1]-y[2];
    tmp[0][1] = x[2]-x[1];
    tmp[1][0] = y[2]-y[0];
    tmp[1][1] = x[0]-x[2];
    tmp[2][0] = y[0]-y[1];
    tmp[2][1] = x[1]-x[0];
    double D = fabs(tmp[2][1]*tmp[1][0] - tmp[1][1]*tmp[2][0]);
    int index = 9*e;
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
        M->idR[index] = s[i];
        M->idC[index] = s[j];
        M->val[index] = (i == j) ? (D/12.) : (D/24.);
        K->idR[index] = s[i];
        K->idC[index] = s[j];
        K->val[index] = 1./(2.*D) * (tmp[i][0]*tmp[j][0] + tmp[i][1]*tmp[j][1]);
        index++;
      }
    }
    
  }
  
  M->nnz = mesh->elements->num * 9;
  K->nnz = mesh->elements->num * 9;
  
  printf("   ... assembling: sorting\n");
  
  sortMatrixSp(M);
  sortMatrixSp(K);
  
  printf("   ... assembling: adding\n");
  
  int indexCurrent = 0;
  
  for(int index=1; index<mesh->elements->num * 9; ++index){
    if((M->idR[index] == M->idR[indexCurrent]) &&
       (M->idC[index] == M->idC[indexCurrent])){
      M->val[indexCurrent] += M->val[index];
      K->val[indexCurrent] += K->val[index];
    }
    else{
      indexCurrent++;
      M->idR[indexCurrent] = M->idR[index];
      M->idC[indexCurrent] = M->idC[index];
      M->val[indexCurrent] = M->val[index];
      K->idR[indexCurrent] = K->idR[index];
      K->idC[indexCurrent] = K->idC[index];
      K->val[indexCurrent] = K->val[index];
    }
  }
  
  M->nnz = indexCurrent+1;
  K->nnz = indexCurrent+1;
  
  // ------------------------------------------------------------------------------
  // 4. Build LINEAR SYSTEM
  // ------------------------------------------------------------------------------
  
  printf("== build linear sytem\n");
  
  // Build matrix "A" (= M+K)
  MatrixSp *A = createMatrixSp(numDof, numDof, M->nnz);
  A->nnz = M->nnz;
  for(int n=0; n<A->nnz; ++n){
    A->idR[n] = M->idR[n];
    A->idC[n] = M->idC[n];
    A->val[n] = M->val[n] + K->val[n];
  }
  
  // Build vector "b" (= M.f with sparce matrix/vector product)
  Vector *b = createVector(numDof);
  for(int i=0; i<numDof; ++i){
    b->val[i] = 0.;
  }
  for(int n=0; n<M->nnz; n++){
    int i = M->idR[n];
    int j = M->idC[n];
    b->val[i] += M->val[n] * f->val[j];
  }
  
  // ------------------------------------------------------------------------------
  // 5. Solve PROBLEM with Jacobi method (https://en.wikipedia.org/wiki/Jacobi_method#Algorithm)
  // ------------------------------------------------------------------------------
  
  printf("== run solver 'jacobi'\n");
  
  // Parameters
  double tol = 1e-6;
  int maxit  = 100000;
  
  // Matrices for Jacobi solver
  Vector   *diagA = createVector(numDof);
  int N_size 	  = A->nnz - numDof;
	printf("N_size: %d, numDof: %d\n", N_size, numDof);
  MatrixSp *N     = createMatrixSp(numDof, numDof, A->nnz - numDof);
  Vector   *Nu    = createVector(numDof);
  for(int n=0; n<A->nnz; n++){
    int i = A->idR[n];
    int j = A->idC[n];
    if(i == j){
      diagA->val[i] = A->val[n];
    }
    else{
      int index = N->nnz;
      N->idR[index] = i;
      N->idC[index] = j;
      N->val[index] = -A->val[n];
      N->nnz = index+1;
    }
  }
  
  // Jacobi solver
  Vector *u = solNum;
  int it = 0;
  double residual  	= 1e9;
	double* diagA_val = diagA->val;
	double* b_val 		= b->val;
	double* Nu_val 		= Nu->val;
	double* u_val 		= u->val;
	double* N_val 		= N->val;
	int* N_idR 				= N->idR;
	int* N_idC 				= N->idC;	
	int nnz 					= N->nnz;

//#pragma acc data pcopy(N_idR[:N_size], N_idC[:N_size], N_val[:N_size], diagA_val[:numDof], b_val[:numDof], Nu_val[:numDof], u_val[:numDof]) 
{
  while (it < maxit && residual > tol){
    printf("residual: %f, tol: %f\n", residual, tol);    
    // Compute N.u
//#pragma acc parallel loop
    for(int i=0; i<numDof; ++i){	//Initialisation des composantes a 0
      Nu_val[i] = 0.;
    }
//#pragma acc serial
    for(int n=0; n<nnz; n++){	//Multiplication avec les valeurs non nulles
      int i = N_idR[n];
      int j = N_idC[n];
      Nu_val[i] += N_val[n] * u_val[j];
    }
    
    // Compute residual
    if((it % 1000) == 0){		//Etape accessoire de visualisation
//#pragma acc wait
      residual = 0;
//#pragma acc parallel loop reduction(+:residual) async
      for(int n=0; n<numDof; n++){
        double tmp = diagA_val[n] * u_val[n] - Nu_val[n] - b_val[n];
        residual += tmp*tmp;
      }
      residual = sqrt(residual);
      printf("    %i %e\n", it, residual);
    }
    
    // Update u
//#pragma acc parallel loop
    for(int i=0; i<numDof; ++i){	//Calcul de A^{-1} * (N * u + b)
      u_val[i] = (Nu_val[i] + b_val[i]) / diagA_val[i];
    }
//#pragma acc parallel async
    it++;
  }
//#pragma acc wait
}	//Fin de la clause data

  printf("   -> final iteration: %i (prescribed max: %i)\n", it, maxit);
  printf("   -> final residual: %e (prescribed tol: %e)\n", residual, tol);
  
  // ------------------------------------------------------------------------------
  // 6. Postprocessing
  // ------------------------------------------------------------------------------
  
  printf("== postprocessing\n");
  
  Vector *solErr = createVector(numDof);
  for(int n=0; n<numDof; n++){
    solErr->val[n] = solNum->val[n] - solRef->val[n];
  }
  
  double L2norm = 0.;
  for(int n=0; n<M->nnz; n++){
    int i = M->idR[n];
    int j = M->idC[n];
    L2norm += solErr->val[i] * M->val[n] * solErr->val[j];
  }
  printf("   -> L2-error = %e\n", sqrt(L2norm));
  
  saveToMsh(solNum->val, mesh, "solNum", "solNum.msh");
  saveToMsh(solErr->val, mesh, "solErr", "solErr.msh");
  
  // ------------------------------------------------------------------------------
  // 7. Free space
  // ------------------------------------------------------------------------------
  
  destroyMesh(mesh);
  destroyMatrixSp(M);
  destroyMatrixSp(K);
  destroyMatrixSp(A);
  destroyMatrixSp(N);
  destroyVector(b);
  destroyVector(solNum);
  destroyVector(solRef);
  destroyVector(solErr);
  destroyVector(f);
  destroyVector(diagA);
  destroyVector(Nu);
  
  return 0;
}
