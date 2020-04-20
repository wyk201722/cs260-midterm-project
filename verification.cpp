#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "verification.h"
#include <math.h>


using namespace std;


int main(int argc, char** argv)
{   
        int n = 4096;
        int b;      
        double *A1,*B1, *C1,*C2;
        double t0, t1,t2,t3;
	double gflops;

        A1 = (double *) malloc(sizeof(double) * n * n);
        B1 = (double *) malloc(sizeof(double) * n * n);
	C1 = (double *) malloc(sizeof(double) * n * n);
	C2 = (double *) malloc(sizeof(double) * n * n);

        if ( randomize_matrix(A1, n, n) ) return -1;
        if ( randomize_matrix(B1, n, n) ) return -1;
	

	for(int i = 0; i < n; i++) {
		for(int j= 0; j < n; j++) {
			C1[i *n + j] = 1;
			C2[i *n + j] = 1;
		}
	}


	t0 = get_sec();
	mm_dac(C1, n, A1, n, B1, n, n);
	t1 = get_sec();
	printf("Elapsed time, Divide and conquer %lf seconds\n", t1 - t0);
	gflops =  pow(pow(2,12),3) / ((t1 - t0) * (pow(10,9)));
	cout<<"Gflops is "<<gflops<<endl;


	t2 = get_sec();
	int i,j,k;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			for( k = 0; k < n; k++) {
				C2[i * n + j] += A1[i * n + k] * B1[k * n + j];
			}
		}
	}
	t3 = get_sec();
	printf("Elapsed time, Linear time %lf seconds\n", t3 - t2);	
	

	if(verify_matrix(C1, C2, n, n) == 0) {
		cout<<"correct"<<endl;	
	} else {
		print_matrix(C1, n ,n);
		print_matrix(C2, n, n);
	}

        free(A1);
        free(B1);
        free(C1);
	free(C2);
	
    return 0;
}



