#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void 		print_matrix		(double* A, int A_l, int A_c);
double* 	transpose		(double* A, int A_l, int A_c);

int main(int argc, char* argv[]){
	int A_l, A_c;
	
	struct timeval t_start, t_stop, t_elapsed;
	float duree;
	
	if(argc == 2){
		A_l = atoi(argv[1]);
		A_c = atoi(argv[1]);
	}else if(argc == 3){
		A_l = atoi(argv[1]);
		A_c = atoi(argv[2]);
	}else{
		A_l = 100;
		A_c = 100;
	}
	
	double* A = malloc(A_l * A_c * sizeof(double));

	for(int i = 0; i < A_l; i++){
		for(int j = 0; j < A_c; j++){
			A[i * A_c + j] = (double) i / (j + 1);
		}
	}
	
	gettimeofday(&t_start, NULL);
	A = transpose(A, A_l, A_c);
        gettimeofday(&t_stop, NULL);
        timersub(&t_stop, &t_start, &t_elapsed);
        duree = t_elapsed.tv_sec + 0.000001 * t_elapsed.tv_usec;
        printf("Temps écoulé:                                   %f\n", duree);
	
	return 0;
}

double* transpose(double* A, int A_l, int A_c){

	int num_blocks = 8;
	int block_l_size = (A_l / num_blocks);
	int block_c_size = (A_c / num_blocks);
	double* A_t = malloc(A_l * A_c * sizeof(double));
	
	for(int i = 0; i < A_l; i += block_l_size){
		for(int j = 0; j < A_c; j += block_c_size){
			for(int k = i; k < i + block_l_size; k++){
				for(int l = j; l < j + block_c_size; l++){
					A_t[l * A_l + k] = A[k * A_c + l];
				}
			}
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
