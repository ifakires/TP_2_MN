#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"

float  mnblas_sasum(const int N, const float *X, const int incX){
    float res = 0.0f;
    for(int i = 0; i < N; i += incX){
        res += X[i];
    }
    return res;

}

double mnblas_dasum(const int N, const double *X, const int incX){
    double res = 0.0f;
    for(int i = 0; i < N; i += incX){
        res += X[i];
    }
    return res;
}

float  mnblas_scasum(const int N, const void *X, const int incX){
    float res = 0.0f;
    for(int i = 0; i < N; i += incX){
        res += ((complexe_float_t*)X) [i].real + ((complexe_float_t*)X) [i].imaginary;;
    }
    return res;
}

double mnblas_dzasum(const int N, const void *X, const int incX){
    double res = 0.0f;
    for(int i = 0; i < N; i += incX){
        res += ((complexe_double_t*)X) [i].real + ((complexe_double_t*)X) [i].imaginary;;
    }
    return res;       
}