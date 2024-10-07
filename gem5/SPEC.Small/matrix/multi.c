#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 2  // Size of the matrices

// Function to multiply two matrices in parallel using OpenMP
void matrixMultiply(double A[N][N], double B[N][N], double C[N][N]) {
    int i, j, k;

    // Set the number of threads
    #pragma omp parallel for private(i, j, k) shared(A, B, C)
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            C[i][j] = 0;
            for (k = 0; k < N; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

int main() {
    double A[N][N], B[N][N], C[N][N];

    // Initialize matrices A and B with random values
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = rand() % 100;
            B[i][j] = rand() % 100;
        }
    }

    // Perform parallel matrix multiplication
    double start_time = omp_get_wtime();  // Start timer
    matrixMultiply(A, B, C);
    double end_time = omp_get_wtime();    // End timer

    printf("Matrix multiplication completed in %f seconds.\n", end_time - start_time);

    return 0;
}
