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

  
    V [0].real = 3 ;
    V [0].imaginary = 3 ;

    V [1].real = 2 ;
    V [1].imaginary = 2 ;

    V [2].real = 1 ;
    V [2].imaginary = 1 ;
  

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

    vector_init_complexe (vec3, 1.0) ;
    
    vector_print_complexe(vec3);
    printf("Le min se trouve a la position : %d",cblas_icamin(VECSIZE,vec3,1));
}