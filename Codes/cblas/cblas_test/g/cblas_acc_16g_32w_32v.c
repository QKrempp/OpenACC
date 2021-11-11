#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

float * 	custom_cblas_t	(float* A, int A_l, int A_c, float* B, int B_l, int B_c, float* C, int C_l, int C_c, float alpha, float beta);

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
	float* restrict A = malloc(m * n * sizeof(float));
	float* restrict B = malloc(n * q * sizeof(float));
	float* restrict C = malloc(m * q * sizeof(float));
#pragma acc data copyout(A[:m*n], B[:n*q], C[:m*q])
{
#pragma acc parallel loop gang collapse(3)
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			for(int k = 0; k < q; k++){
				A[i * n + j] = (float) i / (j + 1);
				B[j * q + k] = (float) j / (k + 1);
				C[i * q + k] = (float) i / (k + 1);
			}
		}
	}
}
        gettimeofday(&t_start, NULL);
        C = custom_cblas_t(A, m, n, B, n, q, C, m, q, 0.5, 0.5);
        gettimeofday(&t_stop, NULL);
        timersub(&t_stop, &t_start, &t_elapsed);
        duree = t_elapsed.tv_sec + 0.000001 * t_elapsed.tv_usec;
        printf("%f\n", duree);
	return 1;
}



float * custom_cblas_t(float* A, int A_l, int A_c, float* B, int B_l, int B_c, float* C, int C_l, int C_c, float alpha, float beta){

	int A_s = A_l * A_c;
	int B_s = B_l * B_c;
	int C_s = C_l * C_c; 

	float* restrict B_t = malloc(B_s * sizeof(float));
#pragma acc data pcopyin(A[0:A_s], B[0:B_s]) pcopy(C[0:C_s]) pcreate(B_t[0:B_s])
{
#pragma acc parallel loop gang num_gangs(512) num_workers(32) 
	for(int i = 0; i < B_c; i++){
#pragma acc loop worker 
		for(int j = 0; j < B_l; j++){
			B_t[j * B_c + i] = B[i * B_l + j];
		}
	}
#pragma acc parallel loop gang num_gangs(16) num_workers(32) vector_length(32)
	for(int i = 0; i < C_l; i++){
#pragma acc loop worker 
		for(int j = 0; j < C_c; j++){
			C[i * C_c + j] *= beta;
			float c = 0;
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


