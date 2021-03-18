#include "mnblas.h"
#include <stdio.h>
#include <math.h>
#include "complexe.h"

float mnblas_snrm2(const int N, const float *X, const int incX)
{
    float res = 0.0f;
    for (int i = 0; i < N; i++)
    {
        res += X[i] * X[i];
    }
    return sqrtf(res);
}

double mnblas_dnrm2(const int N, const double *X, const int incX)
{
    double res = 0.0;
    for (int i = 0; i < N; i++)
    {
        res += X[i] * X[i];
    }
    return sqrtf(res);
}

float mnblas_scnrm2(const int N, const void *X, const int incX)
{
    float res = 0.0f;
    for (int i = 0; i < N; i++)
    {
        res += (((complexe_double_t *)X)[i].real * ((complexe_double_t *)X)[i].real) + (((complexe_double_t *)X)[i].imaginary * ((complexe_double_t *)X)[i].imaginary);
    }
    return sqrtf(res);
}

double mnblas_dznrm2(const int N, const void *X, const int incX)
{
    double res = 0.0;
    for (int i = 0; i < N; i++)
    {
        res += (((complexe_double_t *)X)[i].real * ((complexe_double_t *)X)[i].real) + (((complexe_double_t *)X)[i].imaginary * ((complexe_double_t *)X)[i].imaginary);
    }
    return sqrtf(res);    
}