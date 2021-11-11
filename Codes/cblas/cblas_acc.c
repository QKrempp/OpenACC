// pgcc -fast -acc -Minfo=all cblas_acc.c -o cblas_acc.pgcc
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

double * 	custom_cblas_t	(double* A, int A_l, int A_c, double* B, int B_l, int B_c, double* C, int C_l, int C_c, double alpha, double beta);
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
        C = custom_cblas_t(A, m, n, B, n, q, C, m, q, 0.5, 0.5);
        gettimeofday(&t_stop, NULL);
        timersub(&t_stop, &t_start, &t_elapsed);
        duree = t_elapsed.tv_sec + 0.000001 * t_elapsed.tv_usec;
        printf("Temps écoulé:                                   %f\n", duree);
        if(argc == 2){
                int m_ord = floor(log(m)/log(10));
                int n_ord = floor(log(n)/log(10));
                int q_ord = floor(log(q)/log(10));
                float m_rad = m / pow(10, m_ord);
                float n_rad = n / pow(10, n_ord);
                float q_rad = n / pow(10, q_ord);
                int ops_ord_p = floor(log(m_rad * q_rad * (3 * n + 2))/log(10));
                int ops_ord = ops_ord_p + m_ord + q_ord;
                float ops_rad = m_rad * q_rad * (3 * n + 2) / pow(10, ops_ord_p);
                int flops_ord = floor(log(ops_rad / duree)/log(10));
                float flops_rad = (ops_rad / duree) / pow(10, flops_ord);
                printf("Opérations à virgule flottante effectuées:      %f e %d\n", ops_rad, ops_ord);
                printf("FLOPS théoriques:                               %f e %d\n", flops_rad, flops_ord + ops_ord);
        }
	return 1;
}


double * custom_cblas_t(double* restrict A, int A_l, int A_c, double* restrict B, int B_l, int B_c, double* restrict C, int C_l, int C_c, double alpha, double beta){

	int A_s = A_l * A_c;
	int B_s = B_l * B_c;
	int C_s = C_l * C_c; 

	double* restrict B_t = malloc(B_s * sizeof(double));
#pragma acc data pcopyin(A[0:A_s], B[0:B_s]) pcopy(C[0:C_s]) pcreate(B_t[0:B_s])
{
#pragma acc parallel loop gang num_gangs(512) num_workers(32) 
	for(int i = 0; i < B_c; i++){
#pragma acc loop worker 
		for(int j = 0; j < B_l; j++){
			B_t[j * B_c + i] = B[i * B_l + j];
		}
	}
#pragma acc parallel loop gang num_workers(32) vector_length(32) 
	for(int i = 0; i < C_l; i++){
#pragma acc loop worker 
		for(int j = 0; j < C_c; j++){
			C[i * C_c + j] *= beta;
			double c = 0;
#pragma acc loop vector reduction(+:c)
			for(int k = 0; k < A_c; k++){
				c += alpha * A[i * A_c + k] * B_t[j * B_l + k];
			}
			C[i * C_c + j] += c;
		}
	}	
}
	return C;
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
