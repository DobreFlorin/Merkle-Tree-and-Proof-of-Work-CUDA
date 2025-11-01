/*
 * Tema 2 ASC
 * 2025 Spring
 */
#include "utils.h"
#include <stdlib.h>
#include <string.h>

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double *B, double *x) {
    // Allocate matrices C and D as 2D arrays (array of pointers)
    double **C = (double **)malloc(N * sizeof(double *));
    double **D = (double **)malloc(N * sizeof(double *));
    int i, j, k, iter;

    for (i = 0; i < N; ++i) {
        C[i] = (double *)calloc(N, sizeof(double));
        D[i] = (double *)calloc(N, sizeof(double));
    }

    // Compute C = B * A^T
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            C[i][j] = 0.0;
            for (k = 0; k < N; ++k) {
                C[i][j] += B[i*N + k] * A[j*N + k];
            }
        }
    }

    // Compute D = C^T * A
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            D[i][j] = 0.0;
            for (k = 0; k < N; ++k) {
                D[i][j] += C[k][i] * A[k*N + j];
            }
        }
    }

    // Temporary vectors for the loop
    double *y_loop = (double *)calloc(N, sizeof(double));
    double *x_loop = (double *)calloc(N, sizeof(double));

    for (iter = 0; iter < N; ++iter) {
        // y = C^T * x
        for (i = 0; i < N; ++i) {
            y_loop[i] = 0.0;
            for (k = 0; k < N; ++k) {
                y_loop[i] += C[k][i] * x[k];
            }
        }

        // x = C * y_loop
        for (i = 0; i < N; ++i) {
            x_loop[i] = 0.0;
            for (k = 0; k < N; ++k) {
                x_loop[i] += C[i][k] * y_loop[k];
            }
        }

        // Update x
        memcpy(x, x_loop, N * sizeof(double));
    }

    // Allocate output vector y
    double *y = (double *)calloc(N, sizeof(double));

    // Compute y = D * x
    for (i = 0; i < N; ++i) {
        y[i] = 0.0;
        for (k = 0; k < N; ++k) {
            y[i] += D[i][k] * x[k];
        }
    }

    // Cleanup
    for (i = 0; i < N; ++i) {
        free(C[i]);
        free(D[i]);
    }
    free(C);
    free(D);
    free(y_loop);
    free(x_loop);

    return y;
}
