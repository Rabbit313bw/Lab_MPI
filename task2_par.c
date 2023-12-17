#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>



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


void mul(double* A, double* B, double* C, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
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


void Initialized_matrix(double* A, int n, int m, int my_rank)
{
    srand(time(NULL) + my_rank);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            A[i * m + j] = rand() % 20;
}





int main(int argc, char** argv)
{
    int n = strtol(argv[1], NULL, 10);
    MPI_Init(&argc, &argv);
    int my_rank, comm_sz;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    int q = sqrt(comm_sz);
    int size_block = n / q;
    double* blockA = malloc(size_block * size_block * sizeof(double));
    double* blockB = malloc(size_block * size_block * sizeof(double));
    double* blockC = malloc(size_block * size_block * sizeof(double));
    for (int i = 0; i < size_block; i++)
        for (int j = 0; j < size_block; j++)
            blockC[i * size_block + j] = 0.0;


    double* A = NULL;
    double* B = NULL;
    double* C = NULL;
    Initialized_matrix(blockA, size_block, size_block, my_rank);
    Initialized_matrix(blockB, size_block, size_block, my_rank);
    
    if (my_rank == 0)
    {
        A = malloc(n * n * sizeof(double));
        B = malloc(n * n * sizeof(double));
        C = malloc(n * n * sizeof(double));
    }



    if (my_rank == 0)
    {
        for (int i = 0; i < comm_sz; i++)
        {
            int displ_row = i / q;
            int displ_col = i % q;
            if (i != 0)
            {
                double* buff = malloc(size_block * size_block * sizeof(double));
                MPI_Recv(buff, size_block * size_block, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int j = 0; j < size_block; j++)
                {
                    for (int k = 0; k < size_block; k++)
                    {
                        B[(size_block * displ_row + j) * n + (size_block * displ_col + k)] = buff[j * size_block + k];
                    }
                }
                free(buff);
            }
            else
            {
                for (int j = 0; j < size_block; j++)
                {
                    for (int k = 0; k < size_block; k++)
                    {
                        B[(size_block * displ_row + j) * n + (size_block * displ_col + k)] = blockB[j * size_block + k];
                    }
                }
            }
        }
    }
    else
    {
        MPI_Send(blockB, size_block * size_block, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }



    if (my_rank == 0)
    {
        for (int i = 0; i < comm_sz; i++)
        {
            int displ_row = i / q;
            int displ_col = i % q;
            if (i != 0)
            {
                double* buff = malloc(size_block * size_block * sizeof(double));
                MPI_Recv(buff, size_block * size_block, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int j = 0; j < size_block; j++)
                {
                    for (int k = 0; k < size_block; k++)
                    {
                        A[(size_block * displ_row + j) * n + (size_block * displ_col + k)] = buff[j * size_block + k];
                    }
                }
                free(buff);
            }
            else
            {
                for (int j = 0; j < size_block; j++)
                {
                    for (int k = 0; k < size_block; k++)
                    {
                        A[(size_block * displ_row + j) * n + (size_block * displ_col + k)] = blockA[j * size_block + k];
                    }
                }
            }
        }
    }
    else
    {
        MPI_Send(blockA, size_block * size_block, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }





    double* newblockA = malloc(sizeof(double) * size_block * size_block);
    double* newblockB = malloc(sizeof(double) * size_block * size_block);


    if (my_rank >= q)
    {
        int displ = my_rank / q;
        int send_index = (my_rank - displ) % (displ * q) + q;
        if (send_index < q * displ)
        {
            send_index = my_rank - displ;
        }
        int rec_index = (my_rank + displ) % (displ * q) + q;
        if (rec_index >= q * displ + q || rec_index < q * displ)
        {
            rec_index = my_rank + displ;
        }
        MPI_Sendrecv(blockA, size_block * size_block, MPI_DOUBLE, send_index, 0, newblockA, size_block * size_block, MPI_DOUBLE, rec_index, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // MPI_Send(blockA, size_block * size_block, MPI_DOUBLE, send_index, 0, MPI_COMM_WORLD);
        // MPI_Recv(newblockA, size_block * size_block, MPI_DOUBLE, rec_index, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
        for (int i = 0; i < size_block; i++)
            for (int j = 0; j < size_block; j++)
            {
                newblockA[i * size_block + j] = blockA[i * size_block + j];
            }
    }
 

    if (my_rank % q != 0)
    {
        int displ = my_rank % q;
        int send_index = my_rank - displ * q;
        if (send_index <= 0)
        {
            send_index += q * q;
        }
        int rec_index = my_rank + displ * q;
        if (rec_index > comm_sz)
        {
            rec_index -= q * q;
        }
        MPI_Sendrecv(blockB, size_block * size_block, MPI_DOUBLE, send_index, 2, newblockB, size_block * size_block, MPI_DOUBLE, rec_index, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // MPI_Send(blockB, size_block * size_block, MPI_DOUBLE, send_index, 1, MPI_COMM_WORLD);
        // MPI_Recv(newblockB, size_block * size_block, MPI_DOUBLE, rec_index, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
        for (int i = 0; i < size_block; i++)
            for (int j = 0; j < size_block; j++)
            {
                newblockB[i * size_block + j] = blockB[i * size_block + j];
            }
    }

    double time_start = MPI_Wtime();
    mul(newblockA, newblockB, blockC, size_block);
    double time_finish = MPI_Wtime();
    double local_time = time_finish - time_start;


    int displ = my_rank / q;
    int send_index_A = my_rank - 1;
    if (send_index_A < q * displ)
    {
        send_index_A = my_rank - 1 + q;
    }
    int rec_index_A = my_rank + 1;
    if (my_rank == q * displ + q - 1 || my_rank == q - 1)
    {
        rec_index_A = displ * q;
    }
    displ = my_rank % q;
    int send_index_B = my_rank - q;
    if (send_index_B < 0)
    {
        send_index_B += q * q;
    }
    int rec_index_B = my_rank + q;
    if (rec_index_B >= comm_sz)
    {
        rec_index_B -= q * q;
    }
    


    for (int i = 1; i < q; i++)
    {
        double* tmp = newblockA;
        newblockA = blockA;
        blockA = tmp;
        tmp = newblockB;
        newblockB = blockB;
        blockB = tmp;
        MPI_Sendrecv(blockA, size_block * size_block, MPI_DOUBLE, send_index_A, 1, newblockA, size_block * size_block, MPI_DOUBLE, rec_index_A, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(blockB, size_block * size_block, MPI_DOUBLE, send_index_B, 2, newblockB, size_block * size_block, MPI_DOUBLE, rec_index_B, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        time_start = MPI_Wtime();
        mul(newblockA, newblockB, blockC, size_block);
        time_finish = MPI_Wtime();
        local_time += time_finish - time_start;
    }



    double time_elapsed;
    MPI_Reduce(&local_time, &time_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);




    if (my_rank == 0)
    {
        for (int i = 0; i < comm_sz; i++)
        {
            int displ_row = i / q;
            int displ_col = i % q;
            if (i != 0)
            {
                double* buff = malloc(size_block * size_block * sizeof(double));
                MPI_Recv(buff, size_block * size_block, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int j = 0; j < size_block; j++)
                {
                    for (int k = 0; k < size_block; k++)
                    {
                        C[(size_block * displ_row + j) * n + (size_block * displ_col + k)] = buff[j * size_block + k];
                    }
                }
                free(buff);
            }
            else
            {
                for (int j = 0; j < size_block; j++)
                {
                    for (int k = 0; k < size_block; k++)
                    {
                        C[(size_block * displ_row + j) * n + (size_block * displ_col + k)] = blockC[j * size_block + k];
                    }
                }
            }
        }
    }
    else
    {
        MPI_Send(blockC, size_block * size_block, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }




    if (my_rank == 0)
    {
        printf("Time is %lf\n", time_elapsed);
    }


    if (my_rank == 0)
    {   
        printf("Matrix A is \n");
        Print_matrix(A, n, n);
        printf("Matrix B is \n");
        Print_matrix(B, n, n);
        double* seq_res = malloc(n * n * sizeof(double));
        seq_mul(A, B, seq_res, n);
        printf("My result is \n");
        Print_matrix(C, n, n);
        printf("Correct result is \n");
        Print_matrix(seq_res, n, n);
        free(seq_res);
    }


    





    free(A);
    free(B);
    free(blockA);
    free(blockB);
    free(newblockA);
    free(newblockB);
    free(C);
    MPI_Finalize();
    return 0;
}