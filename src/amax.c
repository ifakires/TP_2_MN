#include "mnblas.h"
#include <stdio.h>
#include "complexe.h"

float abs_f2(float x){
    if(x < 0){
        return -x;
    }
    else{return x;}
}
double abs_d2(double x){
    if(x < 0){
        return -x;
    }
    else{return x;}
}

float asum_absf2(complexe_float_t x){
    return abs_f2(x.real)+abs_f2(x.imaginary);
}
double asum_absd2(complexe_double_t x){
    return abs_d2(x.real)+abs_d2(x.imaginary);
}





int cblas_isamax (const int n, const float *x, const int incx)
{
    if(incx < 0 || n < 0){
        return 0;
    }

    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
    float max = -1;

    for(i = 0; i < n; i += incx){
        /*
        if(x[i] == NULL){
            return i;
        }
        */
        if(x[i] > max){
            max = x[i];
        }
    }

    for(j = 0; j < n; j += incx){
        if(x[j] == max){
           return j;
        }
    }
    return -1;
}

int cblas_idamax (const int n, const double *x, const int incx)
{
    if(incx < 0 || n < 0){
        return 0;
    }

    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
    double max = -1;

    for(i = 0; i < n; i += incx){
        /*
        if(x[i] == NULL){
            return i;
        }
        */
        if(x[i] > max){
            max = x[i];
        }
    }

    for(j = 0; j < n; j += incx){
        if(x[j] == max){
           return j;
        }
    }
    return -1;
}

int cblas_icamax (const int n, const void *x, const int incx)
{
     if(incx < 0 || n < 0){
        return 0;
    }

    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
    float max = -1;

    for(i = 0; i < n; i += incx){
        /*
        if(((complexe_float_t*)x) [i] == NULL){
            return i;
        }
        */
        if(asum_absf2(((complexe_float_t*)x) [i]) > max){
            max = asum_absf2(((complexe_float_t*)x) [i]);
        }
    }

    for(j = 0; j < n; j += incx){
        if(asum_absf2(((complexe_float_t*)x) [j]) == max){
           return j;
        }
    }
    return -1;
}


int cblas_izamax (const int n, const void *x, const int incx)
{
    if(incx < 0 || n < 0){
        return 0;
    }

    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
    double max = -1;

    for(i = 0; i < n; i += incx){
        /*
        if(((complexe_double_t*)x) [i] == NULL){
            return i;
        }
        */
        if(asum_absd2( ((complexe_double_t*)x) [i]) > max){
            max = asum_absd2(((complexe_double_t*)x) [i]);
        }
    }

    for(j = 0; j < n; j += incx){
        if(asum_absd2( ((complexe_double_t*)x) [i] ) == max){
           return j;
        }
    }
    return -1;
}