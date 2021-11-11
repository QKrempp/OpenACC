#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]){
	int x = 0;
	int* arr = malloc(10 * sizeof(int));
	for(int i = 0; i < 10; i++){
		arr[i] = i;
	}
#pragma acc parallel copyin(arr[:10]) copyout(x)
{
#pragma acc loop reduction(+:x)
	for(int i = 0; i < 10; i++){
		x += arr[i];
	}
}	
	printf("%d\n", x);
	return 0;
}
