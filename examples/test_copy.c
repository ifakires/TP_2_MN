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

    vector_init_complexe (vec3, 1.0) ;
    vector_init_complexe (vec4, 2.0) ;

    printf("Affichage des deux vecteurs : \n");
    vector_print_complexe(vec3);
    vector_print_complexe(vec4);
    mncblas_ccopy(VECSIZE,vec3,1,vec4,1);
    printf("Affichage des deux vecteurs apres le swap : \n");
    vector_print_complexe(vec3);
    vector_print_complexe(vec4);


}
