#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"

void cblas_saxpy (const int n, const float a, const float *x, const int incx, float *y, const int incy)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    for (i = 0 ; i < n ; i += incx)
    {
      y[j] += a*x [i] + y [j] ;
      j+=incy ;
    }
    return;
}


void cblas_daxpy (const int n, const double a, const double *x, const int incx, double *y, const int incy)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    for (i = 0 ; i < n ; i += incx)
    {
      y[j] += a*x [i] + y [j] ;
      j+=incy ;
    }
}

void cblas_caxpy (const int n, const void *a, const void *x, const int incx, void *y, const int incy)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    for (i = 0 ; i < n ; i += incx)
    {
      ((complexe_float_t*)y) [j] =  add_complexe_float( mult_complexe_float((*(complexe_float_t*)a),((complexe_float_t*)x) [i]), ((complexe_float_t*)y) [j]) ;
      j+=incy ;
    }
    return;
}

void cblas_zaxpy (const int n, const void *a, const void *x, const int incx, void *y, const int incy)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    for (i = 0 ; i < n ; i += incx)
    {
      ((complexe_double_t*)y) [j] = add_complexe_double( mult_complexe_double( (*(complexe_double_t*)a) ,((complexe_double_t*)x) [i]), ((complexe_double_t*)y) [j]) ;
      j+=incy ;
    }
    return;
}