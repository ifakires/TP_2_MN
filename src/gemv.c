#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"

void mncblas_sgemv(const MNCBLAS_LAYOUT layout,
                   const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const float alpha, const float *A, const int lda,
                   const float *X, const int incX, const float beta,
                   float *Y, const int incY)
{
    int i, j;
    register float res;

    for (i = 0; i != M; i++)
    {
        res = 0.0f;
        res = Y[i];
        res = Y[i] * beta;
        for (j = 0; j != N; j++)
            res += A[i + j * lda] * X[j];
        Y[i] = res;
    }
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout,
                   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const double alpha, const double *A, const int lda,
                   const double *X, const int incX, const double beta,
                   double *Y, const int incY)
{
    int i, j;
    register double res;

    for (i = 0; i != M; i++)
    {
        res = 0.0;
        res = Y[i];
        res = Y[i] * beta;
        for (j = 0; j != N; j++)
            res += A[i + j * lda] * X[j];
        Y[i] = res;
    }
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout,
                   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const void *alpha, const void *A, const int lda,
                   const void *X, const int incX, const void *beta,
                   void *Y, const int incY)
{
    int i, j;
    register complexe_float_t res;
    res.real = 0.0;
    res.imaginary = 0.0;
    res = ((complexe_float_t *)Y)[i];
    res = mult_complexe_float(((complexe_float_t *)Y)[i], *(complexe_float_t *)beta);
    for (j = 0; j != N; j++)
    {
        res = add_complexe_float(res, mult_complexe_float(((complexe_float_t *)X)[j], ((complexe_float_t *)A)[i * j + lda]));
    }
    ((complexe_float_t *)Y)[i] = res;
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout,
                   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const void *alpha, const void *A, const int lda,
                   const void *X, const int incX, const void *beta,
                   void *Y, const int incY)
{
    int i, j;
    register complexe_double_t res;
    res.real = 0.0;
    res.imaginary = 0.0;
    res = ((complexe_double_t *)Y)[i];
    res = mult_complexe_double(((complexe_double_t *)Y)[i], *(complexe_double_t *)beta);
    for (j = 0; j != N; j++)
    {
        res = add_complexe_double(res, mult_complexe_double(((complexe_double_t *)X)[j], ((complexe_double_t *)A)[i * j + lda]));
    }
    ((complexe_double_t *)Y)[i] = res;
}