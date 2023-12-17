#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"



void Print_matrix(double* A, int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
            printf("%.3lf ", A[i * m + j]);
        printf("\n");
    }
}


void matrix_mul(double* A, double* vec, double* res, int start, int finish)
{

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
            A[i * m + j] = (double)(rand() % 20);
}

void Initialized_vector(double* v, int n)
{
    for (int i = 0; i < n; i++)
        v[i] = (double)(rand() % 20);
}

void Print_vector(double* v, int n)
{
    for (int i = 0; i < n; i++)
        printf("%.3lf ", v[i]);
    printf("\n");
}

int main(int argC, char *argV[])
{
    int n, m;
    double* A;
    double* vec;
    int my_rank;
    int comm_sz;
    MPI_Init(&argC, &argV);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


    n = strtol(argV[1], NULL, 10);
    m = strtol(argV[2], NULL, 10);
    A = malloc(n * m * sizeof(double));
    vec = malloc(m * sizeof(double));

    Initialized_matrix(A, n, m);
    Initialized_vector(vec, m);


    int local_cols = m / comm_sz;
    int* send_counts = malloc(comm_sz * sizeof(int));
    int* displs = malloc(comm_sz * sizeof(int));
    int offset = 0;

    for (int i = 0; i < comm_sz; i++)
    {
        if (i == comm_sz - 1)
            send_counts[i] = m - local_cols * (comm_sz - 1);
        else
            send_counts[i] = local_cols;
        displs[i] = offset;
        offset += local_cols;
    }


    double* res = NULL;
    if (my_rank == 0)
        res = malloc(n * sizeof(double));


    double local_time_start, local_time_finish;
    double local_time = 0.0;
    
    for(int i = 0; i < n; i++)
    {
        double local_sum = 0.0;
        local_time_start = MPI_Wtime();
        for (int j = 0; j < send_counts[my_rank]; j++)
        {
            local_sum += vec[j + displs[my_rank]] * A[i * m + j + displs[my_rank]];
        }
        local_time_finish = MPI_Wtime();
        local_time += local_time_finish - local_time_start;
        MPI_Reduce(&local_sum, &res[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    
    double local_time_elapsed = local_time;
    double time_elapsed;
    MPI_Reduce(&local_time_elapsed, &time_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        printf("Time is %lf\n", time_elapsed);
    }

    // if (my_rank == 0)
    // {
    //     printf("My result is \n");
    //     Print_vector(res, n);
    //     double* res_seq = malloc(n * sizeof(double));
    //     seq_mul(A, vec, res_seq, n, m);
    //     printf("Correct result is \n");
    //     Print_vector(res_seq, n);
    //     free(res_seq);
    // }
    
    
    free(res);
    free(A);
    free(vec);
    free(send_counts);
    free(displs);


    MPI_Finalize();
    return 0;
}