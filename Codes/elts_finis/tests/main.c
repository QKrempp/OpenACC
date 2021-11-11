// icc -std=c99 matrix.c mesh.c main.c -o myfem

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "matrix.h"
#include "mesh.h"
#include "flt.h"

int main(int argc, char* argv[])
{
  

	struct timeval t_start, t_stop, t_elapsed;
	float duree;

#pragma acc init

  // ------------------------------------------------------------------------------
  // 1. Build MESH
  // ------------------------------------------------------------------------------
  
  // printf("== read mesh\n");
  
  Mesh *mesh = readMsh("square.msh");
  
  // ------------------------------------------------------------------------------
  // 2. Build FIELDS
  // ------------------------------------------------------------------------------
  
  // printf("== build fields\n");
  
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
  
  // printf("== build mass/stiffness matrices\n");
  
  MatrixSp *M = createMatrixSp(numDof, numDof, mesh->elements->num * 9);
  MatrixSp *K = createMatrixSp(numDof, numDof, mesh->elements->num * 9);
  
  // printf("   ... build local matrices\n");
  
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
  
  // printf("   ... assembling: sorting\n");
  
  sortMatrixSp(M);
  sortMatrixSp(K);
  
  // printf("   ... assembling: adding\n");
  
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
  
  // printf("== build linear sytem\n");
  
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
  

  // printf("== run solver 'jacobi'\n");
	
 
  // Parameters
  double tol = 1e-3;
  int maxit  = 100000;
  
  // Matrices for Jacobi solver
  Vector   *diagA = createVector(numDof);
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
  double residual = 1e9;

	// printf("Renommage des tableaux\n");

	double* restrict Nu_val 		= Nu->val;
	double* restrict N_val 			= N->val;
	double* restrict u_val 			= u->val;
	double* restrict diagA_val 	= diagA->val;
	double* restrict b_val 			= b->val;
	int* restrict N_idR 				= N->idR;
	int* restrict N_idC 				= N->idC;
	int N_nnz 									= N->nnz;

	
	int* restrict N2_nbvR				= calloc(numDof, sizeof(int));

	for(int n = 0; n < N_nnz; n++){
		N2_nbvR[N_idR[n]]++;
	}

	
	int* restrict N2_idzR 			= malloc((numDof + 1) * sizeof(int));
	int* tmp										= calloc(numDof, sizeof(int));

	N2_idzR[0] = 0;
	for(int n = 0; n < numDof; n++){
		N2_idzR[n + 1] = N2_idzR[n] + N2_nbvR[n];
	}
	free(N2_nbvR);

	double* restrict N2_val 		= malloc(N_nnz * sizeof(double));
	int* restrict N2_idC 				= malloc(N_nnz * sizeof(int));


	for(int n = 0; n < N_nnz; n++){
		int ind = tmp[N_idR[n]] + N2_idzR[N_idR[n]];
		N2_val[ind] = N_val[n];
		N2_idC[ind] = N_idC[n];
		tmp[N_idR[n]]++;
	}
	free(tmp);

	// printf("Arrivée dans la zone OpenACC\n");

	gettimeofday(&t_start, NULL); 

#pragma acc data copyin(N2_val[:N_nnz], diagA_val[:numDof], b_val[:numDof], N2_idzR[:numDof + 1], N2_idC[:N_nnz]) copy(u_val[:numDof]) create(Nu_val[:numDof])
{
  while (it < maxit && residual > tol){
    
    // Compute N.u
#pragma acc parallel loop async
    for(int i=0; i<numDof; ++i){
      Nu_val[i] = 0.;
    }

#pragma acc parallel loop independent async
		for(int i = 0; i < numDof; i++){
			for(int j = N2_idzR[i]; j < N2_idzR[i + 1]; j ++){
				Nu_val[i] += N2_val[j] * u_val[N2_idC[j]];
			}
		}

//     for(int n=0; n<N_nnz; n++){
//       int i = N_idR[n];
//       int j = N_idC[n];
//       Nu_val[i] += N_val[n] * u_val[j];
//     }
    
    // Compute residual
    if((it % 1000) == 0){
      residual = 0;
#pragma acc parallel loop copy(residual) reduction(+:residual) async
      for(int n=0; n<numDof; n++){
        double tmp = diagA_val[n] * u_val[n] - Nu_val[n] - b_val[n];
        residual += tmp*tmp;
      }
#pragma acc wait
      residual = sqrt(residual);
      // printf("    %i %e\n", it, residual);
    }
    
    // Update u
#pragma acc parallel loop async
    for(int i=0; i<numDof; ++i){
      u_val[i] = (Nu_val[i] + b_val[i]) / diagA_val[i];
    }
    it++;
  }
#pragma acc wait
 }

	gettimeofday(&t_stop, NULL);
	timersub(&t_stop, &t_start, &t_elapsed);
	duree = t_elapsed.tv_sec + 0.000001 * t_elapsed.tv_usec;

	free(N2_val);
	free(N2_idzR);
	free(N2_idC);

  // printf("   -> final iteration: %i (prescribed max: %i)\n", it, maxit);
  // printf("   -> final residual: %e (prescribed tol: %e)\n", residual, tol);
 
		
	printf("%f\n", duree);
	// printf("Flops arithmétiques:	");
	// prtFlt(divFlt(mltFlt(crtFlt((double)(it * numDof)), crtFlt((double)(3 * numDof + 1))), crtFlt((double)duree)));
	// printf("\n");
	// printf("Flops informatiques:	");
	// prtFlt(divFlt(crtFlt((double)(it * (4 * numDof + 2 * N_nnz) + (it / 1000) * (5 * numDof))), crtFlt((double)duree)));
	// printf("\n");
 
  // ------------------------------------------------------------------------------
  // 6. Postprocessing
  // ------------------------------------------------------------------------------
  
  // printf("== postprocessing\n");
  
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
  // printf("   -> L2-error = %e\n", sqrt(L2norm));
  
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
