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
    A = NULL;
    if (my_rank == 0) 
    {
        A = malloc(n * m * sizeof(double));
        Initialized_matrix(A, n, m);
    }
    vec = malloc(m * sizeof(double));


    Initialized_vector(vec, m);


    int local_rows = n / comm_sz;
    int* send_counts = malloc(comm_sz * sizeof(int));
    int* displs = malloc(comm_sz * sizeof(int));
    int* displs_matrix = malloc(comm_sz * sizeof(int));
    int* send_counts_matrix = malloc(comm_sz * sizeof(int));
    int offset = 0;

    for (int i = 0; i < comm_sz; i++)
    {
        if (i == comm_sz - 1)
            send_counts[i] = n - local_rows * (comm_sz - 1);
        else
            send_counts[i] = local_rows;
        send_counts_matrix[i] = send_counts[i] * m;
        displs[i] = offset;
        displs_matrix[i] = displs[i] * m;
        offset += local_rows;
    }

    double* local_matrix = malloc(send_counts[my_rank] * m * sizeof(double));

    MPI_Scatterv(A, send_counts_matrix, displs_matrix, MPI_DOUBLE, local_matrix, send_counts_matrix[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double local_time_start = MPI_Wtime();

    double* local_res = malloc(send_counts[my_rank] * sizeof(double));
    
    for(int i = 0; i < send_counts[my_rank]; i++)
    {
        local_res[i] = 0.0;
        for (int j = 0; j < m; j++)
            local_res[i] += vec[j] * local_matrix[i * m + j];
    }
    double* res = NULL;
    if (my_rank == 0)
        res = malloc(n * sizeof(double));
    

    double local_time_finish = MPI_Wtime();

    MPI_Gatherv(local_res, send_counts[my_rank], MPI_DOUBLE, res, send_counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    double local_time_elapsed = local_time_finish - local_time_start;
    double time_elapsed;
    MPI_Reduce(&local_time_elapsed, &time_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(my_rank == 0)
    {
        printf("Time is %lf\n", local_time_elapsed);
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
    free(local_matrix);
    free(A);
    free(local_res);
    free(vec);
    free(send_counts);
    free(displs);


    MPI_Finalize();
    // free(res);
    return 0;
}