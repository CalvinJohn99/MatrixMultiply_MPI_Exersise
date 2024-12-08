#include <iostream>
#include <mpi.h>

#define TYPENAME int

void matrixMultiply(int N, const TYPENAME* A, const TYPENAME* B, TYPENAME* C, int my_rank, int world_size) {
    int rows_per_process = N / world_size;
    int start_row = my_rank * rows_per_process;
    int end_row = (my_rank + 1) * rows_per_process;

    for (int i = start_row; i < end_row; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                C[i * N + j] += A[i * N + k] * B[k * N + j];
            }
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    int N;
    if (rank == 0) {
        std::cin >> N;  // Read the size of the matrix from stdin
    }
    
    // Broadcast matrix size to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    TYPENAME* A = new TYPENAME[N * N];
    TYPENAME* B = new TYPENAME[N * N];
    TYPENAME* C = new TYPENAME[N * N]{}; // Result matrix

    // Only process 0 reads the matrices A and B
    if (rank == 0) {
        // Read matrix A
        for (int i = 0; i < N * N; i++) {
            std::cin >> A[i];
        }
        // Read matrix B
        for (int i = 0; i < N * N; i++) {
            std::cin >> B[i];
        }
    }

    // Broadcast matrices A and B to all processes
    MPI_Bcast(A, N * N, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(B, N * N, MPI_INT, 0, MPI_COMM_WORLD);

    // Each process computes its portion of C
    matrixMultiply(N, A, B, C, rank, size);

    if (rank == 0) {
        // Process 0 collects the results from other processes

        for (int neighbour = 1; neighbour < size; neighbour++) {
            int rows_per_process = N / size;
            int start_row = neighbour * rows_per_process;
            int end_row = (neighbour + 1) * rows_per_process;

            MPI_Recv(&C[start_row * N], (end_row - start_row) * N, MPI_INT, neighbour, 100 + neighbour, MPI_COMM_WORLD, &status);
        }

    } else {
        // Send the computed part of C to process 0
        int rows_per_process = N / size;
        int start_row = rank * rows_per_process;
        int end_row = (rank + 1) * rows_per_process;

        MPI_Send(&C[start_row * N], (end_row - start_row) * N, MPI_INT, 0, 100 + rank, MPI_COMM_WORLD);
    }

    // Gather results back to process 0 (optional)
    if (rank == 0) {
        // Print result matrix C
/*         std::cout << "Matrix C:" << std::endl;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                std::cout << C[i * N + j] << " ";
            }
            std::cout << std::endl;
        } */
    }

    // Clean up
    delete[] A;
    delete[] B;
    delete[] C;

    MPI_Finalize();
    return 0;
}






/*
#include "matrixMultiplyMPI.h"

void matrixMultiply_MPI(int N, const floatType* A, const floatType* B, floatType* C, int* flags, int flagCount)
{
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Status status;
    // printf("Hello from process %d out of %d\n", my_rank, world_size);

    int ROWSperrank = int(float(N) / float(world_size));
    int myfirstROW = my_rank*ROWSperrank;
    int mylastROW = (my_rank + 1) * ROWSperrank;

    memset(C, 0, sizeof(floatType) * N * N);
    int blockSize = 32;

    // #pragma omp parallel for collapse(2) 
    // for (int i = myfirstROW; i < mylastROW; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         for (int k = 0; k < N; k++)
    //         {
    //             //C-order 
    //             C[j * N + i] += A[k * N + i].real() * B[j * N + k].real();
    //         }
    //     }
    // }

    #pragma omp parallel for collapse(2) 
    for (int j = 0; j < N; j += blockSize)
    {
        for (int i = myfirstROW; i < mylastROW; i += blockSize)
        {
            for (int k = 0; k < N; k += blockSize)
            {
                for (int di = 0; di < blockSize; di++)
                {
                    int I = i + di;
                    for (int dj = 0; dj < blockSize; dj++) 
                    {
                        int J = j + dj;
                        for (int dk = 0; dk < blockSize; dk++) 
                        {
                            int K = k + dk;
                            C[J * N + I] += A[K * N + I] * B[J * N + K];
                        }
                    }
                }

                //C-order 
                // C[j * N + i] += A[k * N + i].real() * B[j * N + k].real();
            }
        }
    }

    if (my_rank == 0) {
        floatType* C_part = (floatType*)malloc(sizeof(floatType) * N * N);

        // if (C_part == nullptr) {
        //     std::cerr << "Memory allocation failed\n";
        //     MPI_Abort(MPI_COMM_WORLD, -1);
        // }


        for (int neighbour=1; neighbour < world_size; neighbour++) {
            int from_rank;
            MPI_Recv(C_part, N * N, MPI_DOUBLE_COMPLEX, neighbour, 100+neighbour, MPI_COMM_WORLD, &status);
            // MPI_Recv(&from_rank, 1, MPI_INT, neighbour, 101+neighbour, MPI_COMM_WORLD, &status);
            // printf("Process 0 recieved responses from processes %d out of %d\n", from_rank, world_size);
        }

        #pragma omp parallel for collapse(2)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                C[j * N + i] += C_part[j * N + i];
            }
        }

        // #pragma omp parallel for collapse(2)
        // for (int i = 0; i < N; i += blockSize) {
        //     for (int j = 0; j < N; j += blockSize) {
        //         for (int di = 0; di < blockSize; di++) {
        //             int I = i + di;
        //             for (int dj = 0; dj < blockSize; dj++) {
        //                 int J = j + dj;
        //                 C[J * N + I] += C_part[J * N + I];
        //             }
        //         }
        //     }
        // }

        // #pragma omp parallel for collapse(2)
        // for (int ii = 0; ii < N; ii += blockSize) {
        //     for (int jj = 0; jj < N; jj += blockSize) {
        //         // Process each block
        //         for (int i = ii; i < std::min(ii + blockSize, N); i++) {
        //             for (int j = jj; j < std::min(jj + blockSize, N); j++) {
        //                 C[j * N + i] += C_part[j * N + i];
        //             }
        //         }
        //     }
        // }

        free(C_part);

    } else {
        MPI_Send(C, N * N, MPI_DOUBLE_COMPLEX, 0, 100+my_rank, MPI_COMM_WORLD);
        // MPI_Send(&my_rank, 1, MPI_INT, 0, 101+my_rank, MPI_COMM_WORLD);
    }

    MPI_Bcast(C, N * N, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
}
*/