#include <sys/time.h>
#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

using namespace std;


void mm_base(double* c, int n_C,double* a,int n_A,double* b,int n_B,int n) {
   
    /*if(n % 3 == 0) {
    	for (int i = 0; i < n; i += 3)
    	{
        for (int j = 0; j < n; j += 3)
        {
            register double c00 = c[i*n_C + j];
            register double c01 = c[i*n_C + j+1];
            register double c02 = c[i*n_C + j+2];
            register double c10 = c[(i+1)*n_C + j];
            register double c11 = c[(i+1)*n_C + j+1];
            register double c12 = c[(i+1)*n_C + j+2];
            register double c20 = c[(i+2)*n_C + j];
            register double c21 = c[(i+2)*n_C + j+1];
            register double c22 = c[(i+2)*n_C + j+2];

            for (int k = 0; k < n; k += 3)
            {
                for (int l = 0; l < 3; l++)
                {
                    register double a0 = a[i*n_A + k+l];
                    register double a1 = a[(i+1)*n_A + k+l];
                    register double a2 = a[(i+2)*n_A + k+l];
                    register double b0 = b[(k+l)*n_B + j];
                    register double b1 = b[(k+l)*n_B + j+1];
                    register double b2 = b[(k+l)*n_B + j+2];

                    c00 += a0 * b0;
                    c01 += a0 * b1;
                    c02 += a0 * b2;
                    c10 += a1 * b0;
                    c11 += a1 * b1;
                    c12 += a1 * b2;
                    c20 += a2 * b0;
                    c21 += a2 * b1;
                    c22 += a2 * b2;
                }
            }
            c[i*n_C + j] = c00;
            c[i*n_C + j+1] = c01;
            c[i*n_C + j+2] = c02;
            c[(i+1)*n_C + j] = c10;
            c[(i+1)*n_C + j+1] = c11;
            c[(i+1)*n_C + j+2] = c12;
            c[(i+2)*n_C + j] = c20;
            c[(i+2)*n_C + j+1] = c21;
            c[(i+2)*n_C + j+2] = c22;
            }
        }
    } //else {*/
	for(int i = 0; i < n; i += 2) {
		for(int j = 0; j < n; j += 2) {
			register double c00 = c[i * n_C + j];
			register double c01 = c[i * n_C + j + 1];
			register double c10 = c[(i + 1) * n_C + j];
			register double c11 = c[(i + 1) * n_C + j + 1];
			for(int k = 0; k < n; k += 2) {
				for(int l = 0; l < 2; l++) {
					register double a0 = a[i *n_A + k + l];
					register double a1 = a[(i+1) * n_A+ k + l];
					register double b0 = b[(k+l) * n_B + j];
					register double b1 = b[(k+l) * n_B + j + 1];
					
					c00 += a0 * b0;
					c01 += a0 * b1;
					c10 += a1 * b0;
					c11 += a1 * b1;
				}
		 	}
			 c[i * n_C + j] = c00;
			 c[i * n_C + j + 1] = c01;
			 c[(i+1) * n_C + j] = c10;
			 c[(i+1) * n_C + j + 1] = c11;
			}
		}
        /* for (int i = 0; i < n; ++i)
		for (int k = 0; k < n; ++k)
			for (int j = 0; j < n; ++j)
				c[i*n_C+j] += a[i*n_A+k] * b[k*n_B+j];
	
    	}*/
    
}
void mm_dac(double *C, int n_C,double *A, int n_A,double *B, int n_B,int n,int THRESHOLD)
{ 
	assert((n & (-n)) == n);
	if (n <= THRESHOLD) {
		mm_base(C,n_C, A,n_A, B,n_B, n);
	}
	else {
	#define X(M,r,c) (M + (r*(n_ ## M) + c)*(n/2))
	cilk_spawn mm_dac(X(C,0,0), n_C, X(A,0,0), n_A, X(B,0,0), n_B, n/2,THRESHOLD);
	cilk_spawn mm_dac(X(C,0,1), n_C, X(A,0,0), n_A, X(B,0,1), n_B, n/2,THRESHOLD);
	cilk_spawn mm_dac(X(C,1,0), n_C, X(A,1,0), n_A, X(B,0,0), n_B, n/2,THRESHOLD);
		   mm_dac(X(C,1,1), n_C, X(A,1,0), n_A, X(B,0,1), n_B, n/2,THRESHOLD);
	cilk_sync;
	cilk_spawn mm_dac(X(C,0,0), n_C, X(A,0,1), n_A, X(B,1,0), n_B, n/2,THRESHOLD);
	cilk_spawn mm_dac(X(C,0,1), n_C, X(A,0,1), n_A, X(B,1,1), n_B, n/2,THRESHOLD);
	cilk_spawn mm_dac(X(C,1,0), n_C, X(A,1,1), n_A, X(B,1,0), n_B, n/2,THRESHOLD);
		   mm_dac(X(C,1,1), n_C, X(A,1,1), n_A, X(B,1,1), n_B, n/2,THRESHOLD);
	cilk_sync;
	}
}


int print_matrix(const double *A, const int m, const int n)
{
	int i;
	printf("[");
	for (i = 0; A + i && i < m * n; i++)
	{

		if ((i + 1) % n == 0)
			printf("%5.2f ", A[i]);
		else
			printf("%5.2f, ", A[i]);
		if ((i + 1) % n == 0)
		{
			if (i + 1 < m * n)
				printf(";\n");
		}
	}
	printf("]\n");

    if (i != m * n) return -1;
    return 0;
}


int randomize_matrix(double *A, const int m, const int n)
{
	srand (time(NULL));

	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; A + i * n + j && j < n; j++)
		{
			A[i * n + j] = (double)(rand() % 100) + 0.01 * (rand() % 100);
			if ( (rand() % 2) == 0 )
			{
				A[i * n + j] *= -1.;
			}
		}
        if (j != n) return -1;
	}
    return 0;
}


double get_sec()
{
	struct timeval time;
	gettimeofday(&time, NULL);
	return (time.tv_sec + 1e-6 * time.tv_usec);
}


int matrix_copy(double *C, const double *D, const int m, const int n)
{
	int i;

	for (i = 0; C + i && D + i && i < m * n; i++)
	{
		C[i] = D[i];
	}
    	if (i != m * n) 
    	{
        	if (!(C + i)) 
        	{
            	printf("invalid memory access at %d on input 1, return -1\n", i);
            	return -1;
        	}
        	if (!(D + i))
        	{
            	printf("invalid memory access at %d on input 2, return -2\n", i);
            	return -2;
     	   	}
   	 }

    	return 0;
}

int verify_matrix(const double *C, const double *D, const int m, const int n)
{
	int i;
	double diff;
	for (i = 0; i < m * n; i++)
	{
		diff = abs(C[i] - D[i]);
        	if (diff > 1e-3) 
		{
			break;
		}
	}

	if (diff > 1e-3) 
    	{
        	return -1;
    	}
    	else
    	{
        return 0;
    	}
    
}
