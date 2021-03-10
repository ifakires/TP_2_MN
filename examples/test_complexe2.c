#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe2.h"

#define    NB_FOIS        4194304

#include "flop.h"

int main (int argc, char **argv)
{
 complexe_float_t c1= {1.0, 2.0} ;
 complexe_float_t c2= {3.0, 6.0} ;
 complexe_float_t c3 ;

 complexe_double_t cd1 ;
 complexe_double_t cd2 ;
 complexe_double_t cd3 ;
 complexe_double_t cd4 ;

 unsigned long long int start, end ;
 int i ;

 init_flop () ;
 
 c3 = add_complexe_float (c1, c2) ;

 printf ("c3.r %f c3.i %f\n", c3.real, c3.imaginary) ;

 cd1 = (complexe_double_t) {20.0, -4.0} ;
 cd2 = (complexe_double_t) {3.0, 2.0} ;

 cd3 = add_complexe_double (cd1, cd2) ;

 printf ("cd3.r %f cd3.i %f\n", cd3.real, cd3.imaginary) ;

 start =_rdtsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     cd3 = add_complexe_double (cd1, cd2) ;
   }

 end = _rdtsc () ;

  printf ("apres boucle cd3.real %f cd3.imaginary %f %lld cycles \n", cd3.real, cd3.imaginary, end-start) ;

  cd4 = div_complexe_double(cd1,cd2);
  printf ("%f / %f", cd4.real,cd4.imaginary) ;
  printf("\n");

  calcul_flop ("calcul complexe ", NB_FOIS*4, end-start) ;
  exit (0) ;
}


