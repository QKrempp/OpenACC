//	pgcc -fast -acc -Minfo=all mult_kernel.c -o mult_kernel.x


#include <openacc.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void print_matrix(int * a, int n){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			printf("%d ", a[i * n + j]);
		}
		printf("\n");
	}
}

int main(int argc, char * argv[]){
	int n = atoi(argv[1]);
	printf("Création de matrices de taille %d\n", n);
	int *__restrict a = malloc(n * n * sizeof(int));
	int *__restrict b = malloc(n * n * sizeof(int));
	int *__restrict c = malloc(n * n * sizeof(int));
	time_t start, end;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			a[i * n + j] = i+j;
			b[i * n + j] = i-j;
		}
	}
	time(&start);
#pragma acc kernels copyin(a[0:n*n], b[0:n*n]) copyout(c[0:n*n])
	{
	for(int i = 0; i < n; i++){
#pragma acc loop independent gang
		for(int j = 0; j < n; j++){ 
			int t = 0;
#pragma acc loop independent worker(32)
			for(int k = 0; k < n; k++){
				t += a[i * n + k] * b[k * n + j];
			}
			c[i * n + j] = t;
		}
	}
	}
	time(&end);
	free(a);
	free(b);
	free(c);
	//print_matrix(c, n);
	printf("done in: %f sec\n", (double)(end - start));
	return 0;
}
