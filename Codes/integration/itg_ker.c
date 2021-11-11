/*
 *	gcc -O3 -lm itg.c -o itg.gcc
 *	pgcc -fast -Minfo=all itg.c -o itg.pgcc
 */

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <accel.h>


double 		itg	(double a, double b, int N);
double 		f	(double x);

int main(int argc, char* argv[]){
	struct timeval t_start, t_stop, t_elapsed;
	float duree;
	double a, b, integ;
	int N;
	if(argc == 2){
		a = 0;
		b = M_PI;
		N = atoi(argv[1]);
	}else if(argc == 4){
		a = atof(argv[1]);
		b = atof(argv[2]);
		N = atoi(argv[3]);
	}else{
		a = 0;
		b = M_PI;
		N = 1000;
	}
#pragma acc init device_type(acc_device_nvidia)
	gettimeofday(&t_start, NULL);
	integ = itg(a, b, N);
	gettimeofday(&t_stop, NULL);
	timersub(&t_stop, &t_start, &t_elapsed);
	duree = t_elapsed.tv_sec + 0.000001 * t_elapsed.tv_usec;
//	printf("Intégrale calculée:				%f\n", integ);
	printf("%f\n", duree);
	return 0;
}

double itg(double a, double b, int N){
	double approx = 0.5 * (f(a) + f(b));
	double h = (b - a) / N;
#pragma acc kernels loop reduction(+:approx)
	for(int i = 0; i < N; i++){
		approx += h/(a + i * h);
	}
	return approx;	
}

#pragma acc routine
double f(double x){
	return 1/x;
}
