#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define GET_TIME(now) { \
   struct timeval t; \
   gettimeofday(&t, NULL); \
   now = t.tv_sec + t.tv_usec/1000000.0; \
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

void seq_mul(double* A, double* vec, double* res, int n, int m)
{
    for (int i = 0; i < n; i++)
        res[i] = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            res[i] += A[i * m + j] * vec[j];
}


void Initialized_matrix(double* A, int n, int m)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            A[i * m + j] = rand() % 20;
}

void Initialized_vector(double* v, int n)
{
    for (int i = 0; i < n; i++)
        v[i] = rand() % 20;
}

void Print_vector(double* v, int n)
{
    for (int i = 0; i < n; i++)
        printf("%.3lf ", v[i]);
    printf("\n");
}


int main(int argC, char *argV[])
{
    int n = strtol(argV[1], NULL, 10);
    int m = strtol(argV[2], NULL, 10);
    double* A = malloc(sizeof(double) * n * m);
    double* vec = malloc(sizeof(double) * m);
    double* res = malloc(sizeof(double) * n);
    Initialized_matrix(A, n, m);
    Initialized_vector(vec, m);
    
    // printf("matrix is \n");
    // Print_matrix(A, n, m);
    // printf("vector is \n");
    // Print_vector(vec, m);
    double start, end;
    GET_TIME(start);
    seq_mul(A, vec, res, n, m);
    GET_TIME(end);
    // printf("res is \n");
    // Print_vector(res, n);
    printf("time is %lf\n", end - start);


    free(A);
    free(vec);
    free(res);
    return 0;
}