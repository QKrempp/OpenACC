#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

double* iteration(int T, int N, double* C, double* Cnew);

int main(int argc, char* argv[]){
	int N, T;
	struct timeval t_start, t_stop, t_elapsed;
	float duree;
	if (argc == 3){
		N = atoi(argv[1]);
		T = atoi(argv[2]);
	}else{
		N = 10;
		T = 10;
	}
	double* C = malloc(N * N * sizeof(double));
	double* Cnew = malloc(N * N * sizeof(double));
	gettimeofday(&t_start, NULL);
	C = iteration(T, N, C, Cnew);
	gettimeofday(&t_stop, NULL);
	timersub(&t_stop, &t_start, &t_elapsed);
	duree = t_elapsed.tv_sec + 0.000001 * t_elapsed.tv_usec;
	printf("Temps écoulé: 					%f\n", duree);
	return 0;
}

double* iteration(int T, int N, double* C, double* Cnew){
	double dt = 1./(double)(T-1);
	double dx = 1./(double)(N-1);
	for(int n = 0; n < T; n++){
		for(int i = 1; i < (N-1); i++){
			for(int j = 1; j < (N-1); j++){
				Cnew[N * i + j] = (1 - 4 * dt / (dx * dx)) * C[N * i + j] + dt / (dx * dx) * (C[N * (i + 1) + j] + C[N * (i - 1) + j] + C[N * i + (j + 1)] + C[N * i + (j - 1)]);
			}
		}
	}
	return Cnew;
}
