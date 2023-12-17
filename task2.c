#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define GET_TIME(now) { \
   struct timeval t; \
   gettimeofday(&t, NULL); \
   now = t.tv_sec + t.tv_usec/1000000.0; \
}


void seq_mul(double* A, double* B, double* C, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            C[i * n + j] = 0.0;
            for (int k = 0; k < n; k++)
            {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
}

void Print_matrix(double* A, int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
            printf("%.3lf ", A[i * m + j]);
        printf("\n");
    }
}


void Initialized_matrix(double* A, int n, int m)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            A[i * m + j] = rand() % 20;
}



int main(int argc, char** argv)
{
    int n, m, l;
    n = strtol(argv[1], NULL, 10);

    double* A = malloc(n * n * sizeof(double));
    double* B = malloc(n * n * sizeof(double));
    double* C = malloc(n * n * sizeof(double));

    Initialized_matrix(A, n, n);
    Initialized_matrix(B, n, n);

    double start_time, finish_time;
    GET_TIME(start_time);
    seq_mul(A, B, C, n);
    GET_TIME(finish_time);
    double time = finish_time - start_time;
    printf("Time is %lf\n", time);


    free(A);
    free(B);
    free(C);

    return 0;
}

