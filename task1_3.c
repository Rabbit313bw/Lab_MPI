#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>


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
    A = malloc(n * m * sizeof(double));
    Initialized_matrix(A, n, m);
    vec = malloc(m * sizeof(double));


    Initialized_vector(vec, m);


    
    int s = sqrt(comm_sz);
    int q = sqrt(comm_sz);

    int local_rows = n / s;
    int local_cols = n / q;


    int* displs = malloc(comm_sz * sizeof(int));
    int* displs_vector = malloc(comm_sz * sizeof(int));
    int* displs_res = malloc(comm_sz * sizeof(int));
    int elemnts_in_row = local_cols;
    displs[0] = 0;
    displs_vector[0] = 0;
    displs_res[0] = 0;
    int block_in_row = 1;
    int vector_point = local_cols;
    int current_row = 0;

    for(int i = 1; i < comm_sz; i++)
    {
        if (vector_point < m)
        {
            displs_vector[i] = vector_point;
            vector_point += local_cols;
        }
        else
        {
            displs_vector[i] = 0;
            vector_point = local_cols;
        }
        if (elemnts_in_row < m)
        {
            displs[i] = displs[i - 1] + local_cols;
            displs_res[i] = current_row;
            elemnts_in_row += local_cols;
        }
        else
        {
            displs[i] = local_rows * m * block_in_row;
            current_row += local_rows;
            displs_res[i] = current_row;
            block_in_row++;
            elemnts_in_row = local_cols;
        }
    }

    double* gather_res = NULL;
    if (my_rank == 0)
        gather_res = malloc(sizeof(double) * local_rows * comm_sz);
    double* local_sums = malloc(sizeof(double) * local_rows);
    double local_time_start = MPI_Wtime();
    for (int i = 0; i < local_rows; i++)
    {
        local_sums[i] = 0.0;
        for (int j = 0; j < local_cols; j++)
        {
            local_sums[i] += A[(i * m) + displs[my_rank] + j] * vec[displs_vector[my_rank] + j];
        }
    }
    double local_time_finish = MPI_Wtime();


    MPI_Gather(local_sums, local_rows, MPI_DOUBLE, gather_res, local_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* res = NULL;
    double local_time = 0.0;
    local_time += local_time_finish - local_time_start;
    local_time_start = MPI_Wtime();
    if(my_rank == 0)
    {
        res = malloc(sizeof(double) * n);
        for (int i = 0; i < n; i++)
            res[i] = 0.0;

        
        for (int i = 0; i < comm_sz; i++)
        {
            for (int j = 0; j < local_rows; j++)
                res[displs_res[i] + j] += gather_res[local_rows * i + j];
        }
    }
    local_time_finish = MPI_Wtime();
    local_time += local_time_finish - local_time_start;

    double time_elapsed;
    MPI_Reduce(&local_time, &time_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);



    // if (my_rank == 1)
    // {
    //     Print_matrix(A, n, m);
    //     printf("my rank is %d\n", my_rank);
    //     for (int i = 0; i < local_rows; i++)
    //     {
    //         for (int j = 0; j < local_cols; j++)
    //             printf("%.1lf ", A[(i * m) + displs[my_rank] + j]);
    //         printf("\n");
    //     }
    // }

    // if(my_rank == 2)
    // {
    //     for (int i = 0; i < local_rows; i++)
    //         printf("%d ", displs_res[my_rank] + i);
    //     printf("\n");
    // }
        



    if(my_rank == 0)
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
    
    
    free(A);
    free(vec);
    free(displs);
    free(displs_vector);
    free(displs_res);
    free(local_sums);
    free(gather_res);
    free(res);

    MPI_Finalize();
    // free(res);
    return 0;
}