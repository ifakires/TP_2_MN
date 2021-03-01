
typedef struct {
  float real ;
  float imaginary ;
} complexe_float_t ;

typedef struct {
  double real ;
  double imaginary ;
} complexe_double_t ;

inline complexe_float_t add_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r ;

  r.real = c1.real + c2.real ;
  r.imaginary = c1.imaginary + c2.imaginary ;
  
  return r ;
}

inline complexe_double_t add_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
{
  complexe_double_t r ;

  r.real = c1.real + c2.real ;
  r.imaginary = c1.imaginary + c2.imaginary ;
  
  return r ;
}



inline complexe_float_t mult_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r ;

  r.real = c1.real*c2.real - c1.imaginary*c2.imaginary;
  r.imaginary = c1.real*c2.imaginary + c1.imaginary*c2.real;
  
  return r ;
}

inline complexe_double_t mult_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
{
  complexe_double_t r ;

  r.real = c1.real*c2.real - c1.imaginary*c2.imaginary;
  r.imaginary = c1.real*c2.imaginary + c1.imaginary*c2.real;
  
  return r ;
}


inline complexe_float_t div_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r;
  complexe_float_t nume;
  complexe_float_t denom;

  complexe_float_t conjugue;
  conjugue.real = c2.real;
  conjugue.imaginary = -c2.imaginary;

  nume = mult_complexe_float(c1,conjugue);
  denom = mult_complexe_float(c2,conjugue);

  r.real = nume.real/denom.real;
  r.imaginary =nume.imaginary/denom.real;

  return r;
}

inline complexe_double_t div_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
{
  complexe_double_t r;
  complexe_double_t nume;
  complexe_double_t denom;

  complexe_double_t conjugue;
  conjugue.real = c2.real;
  conjugue.imaginary = -c2.imaginary;

  nume = mult_complexe_double(c1,conjugue);
  denom = mult_complexe_double(c2,conjugue);

  r.real = nume.real/denom.real;
  r.imaginary =nume.imaginary/denom.real;

  return r;
}




