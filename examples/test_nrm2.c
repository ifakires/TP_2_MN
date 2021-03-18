#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"
#include "flop.h"
#include <math.h>

#define VECSIZE 3

typedef float vfloat[VECSIZE];
typedef complexe_float_t vfloat_complexe[VECSIZE];
vfloat vec1, vec2;
vfloat_complexe vec3, vec4;

void vector_init(vfloat V, float x)
{
    register unsigned int i;

    for (i = 0; i < VECSIZE; i++)
        V[i] = x;

    return;
}

void vector_init_complexe(vfloat_complexe V, float x)
{
    register unsigned int i;

    for (i = 0; i < VECSIZE; i++)
    {
        V[i].real = x;
        V[i].imaginary = x;
    }

    return;
}

void vector_print_complexe(vfloat_complexe V)
{
    register unsigned int i;

    for (i = 0; i < VECSIZE; i++)
        printf("(%f , %f) ", V[i].real, V[i].imaginary);
    printf("\n");

    return;
}

int main()
{

    vector_init(vec1, 3.0);
    vector_init_complexe(vec3, 6.0);

    printf("nrm2 float : \n");
    float res = mnblas_snrm2(VECSIZE, vec1, 1);
    printf("%f\n", res);

    printf("nrm2 complex : \n");
    float res2 = mnblas_snrm2(VECSIZE, vec3, 1);
    printf("%f\n", res2);
}
