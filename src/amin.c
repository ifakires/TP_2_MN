#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"

float abs_f1(float x){
    if(x < 0){
        return -x;
    }
    else{return x;}
}
double abs_d1(double x){
    if(x < 0){
        return -x;
    }
    else{return x;}
}

float sum_absf1(complexe_float_t x){
    return abs_f1(x.real)+abs_f1(x.imaginary);
}
double sum_absd1(complexe_double_t x){
    return abs_d1(x.real)+abs_d1(x.imaginary);
}

int cblas_isamin (const int n, const float *x, const int incx)
{
    if(incx < 0 || n < 0){
        return 0;
    }

    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
    float min = -1;

    for(i = 0; i < n; i += incx){
        /*
        if( *(x+i) == NULL){
            return i;
        }
        */
        if(x[i] < min){
            min = x[i];
        }
    }

    for(j = 0; j < n; j += incx){
        if(x[j] == min){
           return j;
        }
    }
    return -1;
}

int cblas_idamin (const int n, const double *x, const int incx)
{
    if(incx < 0 || n < 0){
        return 0;
    }

    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
    double min = -1;

    for(i = 0; i < n; i += incx){
        /*
        if( *(x+i) == NULL){
            return i;
        }
        */
        if(x[i] < min){
            min = x[i];
        }
    }

    for(j = 0; j < n; j += incx){
        if(x[j] == min){
           return j;
        }
    }
    return -1;
}

int cblas_icamin (const int n, const void *x, const int incx)
{
    if(incx < 0 || n < 0){
        return 0;
    }

    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
    float min = sum_absf1(((complexe_float_t*)x) [0]);

    for(i = 0; i < n; i += incx){
        /*
        if(((complexe_float_t*)x) [i] == NULL){
            return i;
        }
        */
        if(sum_absf1(((complexe_float_t*)x) [i]) < min){
            min = sum_absf1(((complexe_float_t*)x) [i]);
        }
    }

    for(j = 0; j < n; j += incx){
        if(sum_absf1(((complexe_float_t*)x) [j]) == min){
           return j;
        }
    }
    return -1;
}

int cblas_izamin (const int n, const void *x, const int incx)
{
    if(incx < 0 || n < 0){
        return 0;
    }

    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
    double min = -1;

    for(i = 0; i < n; i += incx){
       /*
        if( *(x+i) == NULL){
            return i;
        }
        */
        if(sum_absd1(((complexe_double_t*)x) [i]) < min){
            min = sum_absd1(((complexe_double_t*)x) [i]);
        }
    }

    for(j = 0; j < n; j += incx){
        if(sum_absd1(((complexe_double_t*)x) [j]) == min){
           return j;
        }
    }
    return -1;
}