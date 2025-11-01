/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"
#include <cblas.h>
#include <stdlib.h>
#include <string.h>

/* 
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B, double *x) {
    double* C = (double*)malloc(N * N * sizeof(double));
    double* D = (double*)malloc(N * N * sizeof(double));
    double* y_loop = (double*)malloc(N * sizeof(double));
    double* x_loop = (double*)malloc(N * sizeof(double));
    double* y = (double*)malloc(N * sizeof(double));
    int iter;

    // C = B * A^T
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, B, N, A, N, 0.0, C, N);

    // D = C^T * A
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, N, 1.0, C, N, A, N, 0.0, D, N);

    for (iter = 0; iter < N; ++iter) {
        // y_loop = C^T * x
        cblas_dgemv(CblasRowMajor, CblasTrans, N, N, 1.0, C, N, x, 1, 0.0, y_loop, 1);
        // x_loop = C * y_loop
        cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, C, N, y_loop, 1, 0.0, x_loop, 1);
        memcpy(x, x_loop, N * sizeof(double));
    }

    // y = D * x
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, D, N, x, 1, 0.0, y, 1);

    free(C);
    free(D);
    free(y_loop);
    free(x_loop);

    return y;
}
