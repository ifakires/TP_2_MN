#include "mnblas.h"
#include <stdio.h>

/*
float mncblas_sdot(const int N, const float *X, const int incX, 
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float dot = 0.0 ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      dot = dot + X [i] * Y [j] ;
    }

  return dot ;
}
*/

float mncblas_sdot(const int N, const float *X, const int incX, 
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  float dot = 0.0 ;

  
  for (i = 0 ; i < N ; i += incX)
    {
      dot += X [i] * Y [j] ;
      j+=incY ;
    }

  return dot ;
}

double mncblas_ddot(const int N, const double *X, const int incX, 
                 const double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  double dot = 0 ;

  
  for (i = 0 ; i < N ; i += incX)
    {
      dot += X [i] * Y [j] ;
      j+=incY ;
    }

  return dot ;
}

void   mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_double_t res ;

  
  for (i = 0 ; i < N ; i += incX)
    {
      res = add_complexe_double(res , mult_complexe_double(X[i] , Y [j])) ;
      j+=incY ;
    }
  dotc = res;
  return;
}

void   mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_double_t res ;

  
  for (i = 0 ; i < N ; i += incX)
    {
      complexe_double_t conjugue;
      conjugue.real = X[i].real;
      conjugue.imaginary = -Y[j].imaginary;

      res = add_complexe_double(res , mult_complexe_double(conjugue , Y [j])) ;
      j+=incY ;
    }
  dotc = res;
  return;
}

void   mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_float_t res ;

  
  for (i = 0 ; i < N ; i += incX)
    {
      res = add_complexe_float(res , mult_complexe_float(X[i] , Y [j])) ;
      j+=incY ;
    }
  dotc = res;
  return;
}
  
void   mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_float_t res ;

  
  for (i = 0 ; i < N ; i += incX)
    {
      complexe_float_t conjugue;
      conjugue.real = X[i].real;
      conjugue.imaginary = -Y[j].imaginary;

      res = add_complexe_float(res , mult_complexe_float(conjugue , Y [j])) ;
      j+=incY ;
    }
  dotc = res;
  return;
}




