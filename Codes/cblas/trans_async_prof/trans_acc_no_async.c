#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <accel.h>

void 		print_matrix		(double* A, int A_l, int A_c);
double* 	transpose		(double* A, int A_l, int A_c, int threads);

int main(int argc, char* argv[]){
	int A_l, A_c;

	int threads;
	
	struct timeval t_start, t_stop, t_elapsed;
	float duree;
	
	if(argc == 2){
		A_l = atoi(argv[1]);
		A_c = atoi(argv[1]);
		threads = 8;
	}else if(argc == 3){
		A_l = atoi(argv[1]);
		A_c = atoi(argv[1]);
		threads = atoi(argv[2]);
	}else{
		A_l = 100;
		A_c = 100;
		threads = 4;
	}
	
	double* A = malloc(A_l * A_c * sizeof(double));

	for(int i = 0; i < A_l; i++){
		for(int j = 0; j < A_c; j++){
			A[i * A_c + j] = (double) i / (j + 1);
		}
	}
	
	gettimeofday(&t_start, NULL);
	A = transpose(A, A_l, A_c, threads);
        gettimeofday(&t_stop, NULL);
        timersub(&t_stop, &t_start, &t_elapsed);
        duree = t_elapsed.tv_sec + 0.000001 * t_elapsed.tv_usec;
        printf("Temps écoulé:                                   %f\n", duree);
	
	return 0;
}

double* transpose(double* A, int A_l, int A_c, int threads){

	int num_blocks = threads;
	int block_c_size = (A_c / num_blocks);
	double* A_t = malloc(A_l * A_c * sizeof(double));
	
#pragma acc data pcopyin(A[:A_l * A_c]) pcreate(A_t[:A_c * A_l])
{

	for(int i = 0; i < num_blocks; i++){
#pragma acc parallel loop gang num_gangs(256) num_workers(32) 
		for(int k = 0; k < block_c_size; k++){
#pragma acc loop worker
			for(int l = 0; l < A_l; l++){
				A_t[i * block_c_size * A_c + k * A_l + l] = A[l * A_c + i * block_c_size + k];
			}
		}
#pragma acc update self(A_t[i * block_c_size * A_c:A_l * block_c_size]) 
	}
}
	return A_t;
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
