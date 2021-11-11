
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

float * 	custom_cblas	(float* A, int A_l, int A_c, float* B, int B_l, int B_c, float* C, int C_l, int C_c, float alpha, float beta);
float * 	custom_cblas_t	(float* A, int A_l, int A_c, float* B, int B_l, int B_c, float* C, int C_l, int C_c, float alpha, float beta);
void 		print_matrix	(float* A, int A_l, int A_c);

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
		m = 100;
		n = 100;
		q = 100;
	}
	float* A = malloc(m * n * sizeof(float));
	float* B = malloc(n * q * sizeof(float));
	float* C = malloc(m * q * sizeof(float));
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			for(int k = 0; k < q; k++){
				A[i * n + j] = (float) i / (j + 1);
				B[j * q + k] = (float) j / (k + 1);
				C[i * q + k] = (float) i / (k + 1);
			}
		}
	}
	gettimeofday(&t_start, NULL);
	C = custom_cblas_t(A, m, n, B, n, q, C, m, q, 0.5, 0.5);
	gettimeofday(&t_stop, NULL);
	timersub(&t_stop, &t_start, &t_elapsed);
	duree = t_elapsed.tv_sec + 0.000001 * t_elapsed.tv_usec;
	printf("Temps écoulé: 					%f\n", duree);
	if(argc == 2){
		int m_ord = floor(log(m)/log(10));
		int n_ord = floor(log(n)/log(10));
		int q_ord = floor(log(q)/log(10));
		float m_rad = m / pow(10, m_ord);
		float n_rad = n / pow(10, n_ord);
		float q_rad = n / pow(10, q_ord);
		printf("m = %f e %d\n", m_rad, m_ord);
		int ops_ord_p = floor(log(m_rad * q_rad * (3 * n + 2))/log(10));
		int ops_ord = ops_ord_p + m_ord + q_ord;
		float ops_rad = m_rad * q_rad * (3 * n + 2) / pow(10, ops_ord_p);
		int flops_ord = floor(log(ops_rad / duree)/log(10));
		float flops_rad = (ops_rad / duree) / pow(10, flops_ord);
		printf("Opérations à virgule flottante effectuées: 	%f e %d\n", ops_rad, ops_ord);
		printf("FLOPS théoriques: 				%f e %d\n", flops_rad, flops_ord + ops_ord);
	}
	return 1;
}

float * custom_cblas(float* A, int A_l, int A_c, float* B, int B_l, int B_c, float* C, int C_l, int C_c, float alpha, float beta){
	for(int i = 0; i < C_l; i++){
		for(int j = 0; j < C_c; j++){
			C[i * C_c + j] *= beta;
			for(int k = 0; k < A_c; k++){
				C[i * C_c + j] += alpha * A[i * A_c + k] * B[k * B_c + j];
			}
		}
	}
	return C;
}

float * custom_cblas_t(float* A, int A_l, int A_c, float* B, int B_l, int B_c, float* C, int C_l, int C_c, float alpha, float beta){
	float* B_t = malloc(B_l * B_c * sizeof(float));
	for(int i = 0; i < B_c; i++){
		for(int j = 0; j < B_l; j++){
			B_t[j * B_c + i] = B[i * B_l + j];
		}
	}
	for(int i = 0; i < C_l; i++){
		for(int j = 0; j < C_c; j++){
			C[i * C_c + j] *= beta;
			for(int k = 0; k < A_c; k++){
				C[i * C_c + j] += alpha * A[i * A_c + k] * B_t[j * B_l + k];
			}
		}
	}
	return C;
}
void print_matrix(float* A, int A_l, int A_c){
	for(int i = 0; i < A_l; i++){
		for(int j = 0; j < A_c; j++){
			printf("%f ", A[i * A_c + j]);
		}
		printf("\n");
	}
	printf("\n");
}
