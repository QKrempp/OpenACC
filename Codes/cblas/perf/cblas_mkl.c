// icc -O3 -mkl cblas_mkl.c -o cblas_mkl.icc
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "mkl.h"

void 		print_matrix	(double* A, int A_l, int A_c);

int main(int argc, char* argv[]){
	int m, n, q;
	struct timeval t_start, t_stop, t_elapsed;
	float duree;
		
	if(argc == 4){
		m = atoi(argv[1]);
		n = atoi(argv[2]);
		q = atoi(argv[3]);
	}else if(argc == 2){
		m = atoi(argv[1]);
	 	n = atoi(argv[1]);
		q = atoi(argv[1]);
	}else{
		m = 10;
		n = 10;
		q = 10;
	}
	double* restrict A = malloc(m * n * sizeof(double));
	double* restrict B = malloc(n * q * sizeof(double));
	double* restrict C = malloc(m * q * sizeof(double));
#pragma acc data copyout(A[:m*n], B[:n*q], C[:m*q])
{
#pragma acc parallel loop gang collapse(3)
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			for(int k = 0; k < q; k++){
				A[i * n + j] = (double) i / (j + 1);
				B[j * q + k] = (double) j / (k + 1);
				C[i * q + k] = (double) i / (k + 1);
			}
		}
	}
}
        gettimeofday(&t_start, NULL);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, q, 0.5, A, n, B, n, 0.5, C, n);
        gettimeofday(&t_stop, NULL);
        timersub(&t_stop, &t_start, &t_elapsed);
        duree = t_elapsed.tv_sec + 0.000001 * t_elapsed.tv_usec;
	printf("%f\n", duree);
//        printf("Temps écoulé:                                   %f\n", duree);
//        if(argc == 2){
//                int m_ord = floor(log(m)/log(10));
//                int n_ord = floor(log(n)/log(10));
//                int q_ord = floor(log(q)/log(10));
//                float m_rad = m / pow(10, m_ord);
//                float n_rad = n / pow(10, n_ord);
//                float q_rad = n / pow(10, q_ord);
//                int ops_ord_p = floor(log(m_rad * q_rad * (3 * n + 2))/log(10));
//                int ops_ord = ops_ord_p + m_ord + q_ord;
//                float ops_rad = m_rad * q_rad * (3 * n + 2) / pow(10, ops_ord_p);
//                int flops_ord = floor(log(ops_rad / duree)/log(10));
//                float flops_rad = (ops_rad / duree) / pow(10, flops_ord);
//                printf("Opérations à virgule flottante effectuées:      %f e %d\n", ops_rad, ops_ord);
//                printf("FLOPS théoriques:                               %f e %d\n", flops_rad, flops_ord + ops_ord);
//        }
	return 0;
}


void print_matrix(double* A, int A_l, int A_c){
	for(int i = 0; i < A_l; i++){
		for(int j = 0; j < A_c; j++){
			printf("%f ", A[i * A_c + j]);
		}
		printf("\n");
	}
	printf("\n");
}