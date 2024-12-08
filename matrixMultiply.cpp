#include <iostream>
#include <stdio.h>
#include <chrono>
#include <immintrin.h>

#define ALIGN 64
#define TYPENAME double

template<typename T>
void matrixMultiply(int N, const T* A, const T* B, T* C)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                C[i * N + j] += A[i * N + k] * B[k * N + j];
            }
        }
    }
}

int main(int argc, char** argv)
{
    int N;
    std::cin >> N;  // Read the size of the matrix from stdin

    TYPENAME* A = new TYPENAME[N * N];
    TYPENAME* B = new TYPENAME[N * N];
    TYPENAME* C = new TYPENAME[N * N]{};

    // Read matrix A
    for (int i = 0; i < N * N; i++) {
        std::cin >> A[i];
    }

    // Read matrix B
    for (int i = 0; i < N * N; i++) {
        std::cin >> B[i];
    }

    // Call the matrix multiplication function
    matrixMultiply(N, A, B, C);

    // Print result matrix C
    std::cout << "Matrix C:" << std::endl;
/*     for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << C[i * N + j] << " ";
        }
        std::cout << std::endl;
    } */

    // Clean up
    delete[] A;
    delete[] B;
    delete[] C;

    return 0;
}
