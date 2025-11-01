#include "utils.h"
#include <stdlib.h>
#include <string.h>

double* my_solver(int N, double *A, double *B, double *x) {
    int i, j, k, iter;
    double* C = (double*)calloc(N*N, sizeof(double));
    double* D = (double*)calloc(N*N, sizeof(double));
    double* At = (double*)malloc(N*N * sizeof(double));
    double* Ct = (double*)malloc(N*N * sizeof(double));
    double* y_loop = (double*)malloc(N * sizeof(double));
    double* x_work = (double*)malloc(N * sizeof(double));  // Temporary buffer for x
    double* y = (double*)malloc(N * sizeof(double));

    // Precompute At = A^T (correctly transposed)
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            At[i*N + j] = A[j*N + i];  // At[i][j] = A[j][i]
        }
    }

    // Compute C = B * A^T 
    for (i = 0; i < N; ++i) {
        for (k = 0; k < N; ++k) {
            double B_ik = B[i*N + k];  // B[i][k]
            for (j = 0; j < N; ++j) {
                // C[i][j] += B[i][k] * At[k][j]
                C[i*N + j] += B_ik * At[k*N + j];  
            }
        }
    }

    // Precompute Ct = C^T 
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            Ct[i*N + j] = C[j*N + i];  // Ct[i][j] = C[j][i]
        }
    }

    // Compute D = Ct * A (
    for (i = 0; i < N; ++i) {
        for (k = 0; k < N; ++k) {
            double Ct_ik = Ct[i*N + k];  // Ct[i][k]
            for (j = 0; j < N; ++j) {
                D[i*N + j] += Ct_ik * A[k*N + j];  // A[k][j]
            }
        }
    }

    
    memcpy(x_work, x, N * sizeof(double));

    // Optimized loop with pointer swapping
    for (iter = 0; iter < N; ++iter) {
        // y_loop = Ct * x_work (C^T * x)
        for (i = 0; i < N; ++i) {
            double sum = 0.0;
            for (k = 0; k < N; ++k) {
                sum += Ct[i*N + k] * x_work[k];  // Ct[i][k] * x[k]
            }
            y_loop[i] = sum;
        }

        // x_work = C * y_loop
        for (i = 0; i < N; ++i) {
            double sum = 0.0;
            for (k = 0; k < N; ++k) {
                sum += C[i*N + k] * y_loop[k];  // C[i][k] * y_loop[k]
            }
            x_work[i] = sum;
        }
    }

    // Final y = D * x_work
    for (i = 0; i < N; ++i) {
        double sum = 0.0;
        for (k = 0; k < N; ++k) {
            sum += D[i*N + k] * x_work[k];  // D[i][k] * x_work[k]
        }
        y[i] = sum;
    }

    free(C);
    free(D);
    free(At);
    free(Ct);
    free(y_loop);
    free(x_work);

    return y;
}
