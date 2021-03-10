#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"
#include "flop.h"

#define VECSIZE    3

typedef float vfloat [VECSIZE] ;
typedef complexe_float_t vfloat_complexe [VECSIZE] ;
vfloat vec1, vec2 ;
vfloat_complexe vec3, vec4 ;

void vector_init (vfloat V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_init_complexe (vfloat_complexe V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++){
    V [i].real = x ;
    V [i].imaginary = x ;
  }

  return ;
}

void vector_print_complexe (vfloat_complexe V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("(%f , %f) ", V[i].real, V[i].imaginary) ;
  printf ("\n") ;
  
  return ;
}

int main(){

    vector_init(vec1, 3.0) ;
    vector_init_complexe (vec3, 3.0) ;
    vector_init_complexe (vec4, 2.0) ;

    printf("asum : \n");
    float res = mnblas_sasum(VECSIZE, vec1, 1);
    printf("%f\n", res);
    printf("asum complex : \n");
    float res2 = mnblas_scasum(VECSIZE, vec3, 1);
    double res3 = mnblas_scasum(VECSIZE, vec4, 1);
    printf("%f\n", res2);
    printf("%f\n", res3);


}
