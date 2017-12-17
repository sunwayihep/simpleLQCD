/*******************************************************
 *
 * some useful operation on su(3) matrix, vector 
 * and spinor type
 *
 * ****************************************************/
#ifndef __DATATYPE_H_
#define __DATATYPE_H_
#include "global.h"

typedef struct
{
  Float re, im;
}complex;

typedef struct
{
  complex c1, c2, c3;
}su3_vector;

typedef struct
{
  complex c11, c12, c13;
  complex c21, c22, c23;
  complex c31, c32, c33;
}su3;

typedef struct
{
  su3_vector c1,c2,c3,c4;
}spinor;

// gauge and spinor field
typedef su3 su3_field[NT][NZ][NY][NX][4];
typedef spinor spinor_field[NT][NZ][NY][NX];

// multiplication of real and complex number
static inline complex re_complex_mul(complex z, Float f)
{
 
  complex tmp;
  tmp.re = z.re * f;
  tmp.im = z.im * f;
  
  return tmp;

}

// multiplication of two complex number
static inline complex complex_mul(complex a1, complex a2)
{

  complex tmp;
  tmp.re = a1.re*a2.re - a1.im*a2.im;
  tmp.im = a1.im*a2.re + a1.re*a2.im;
  
  return tmp;

}

// multiplication of real number and su3 vector
static inline void re_vector_mul(su3_vector* result, const Float real, 
                                 const su3_vector* const in)
{
 
  result->c1.re= real*(in->c1.re);
  result->c1.im= real*(in->c1.im);
  result->c2.re= real*(in->c2.re);
  result->c2.im= real*(in->c2.im);
  result->c3.re= real*(in->c3.re);
  result->c3.im= real*(in->c3.im);

}

// multiplication of complex number and su3 vector
static inline void cm_vector_mul(su3_vector* result, const complex cm, 
                                 const su3_vector* const in)
{
 
  result->c1.re= cm.re*(in->c1.re) - cm.im*(in->c1.im);
  result->c1.im= cm.im*(in->c1.re) + cm.re*(in->c1.im);
  result->c2.re= cm.re*(in->c2.re) - cm.im*(in->c2.im);
  result->c2.im= cm.im*(in->c2.re) + cm.re*(in->c2.im);
  result->c3.re= cm.re*(in->c3.re) - cm.im*(in->c3.im);
  result->c3.im= cm.im*(in->c3.re) + cm.re*(in->c3.im);

}

// result = i * in
static inline void I_vector_mul(su3_vector*  result, 
                                const su3_vector* const in)
{

  result->c1.re= -(in->c1.im);
  result->c1.im= (in->c1.re);
  result->c2.re= -(in->c2.im);
  result->c2.im= (in->c2.re);
  result->c3.re= -(in->c3.im);
  result->c3.im= (in->c3.re);

}

// result = -i * in
static inline void mI_vector_mul(su3_vector* result, 
                                 const su3_vector* const in)
{

  result->c1.re= (in->c1.im);
  result->c1.im= -(in->c1.re);
  result->c2.re= (in->c2.im);
  result->c2.im= -(in->c2.re);
  result->c3.re= (in->c3.im);
  result->c3.im= -(in->c3.re);

}

// result += input
static inline void vector_add_assign(su3_vector* result, 
                                     const su3_vector* const input)
{

  result->c1.re += input->c1.re;
  result->c1.im += input->c1.im;
  result->c2.re += input->c2.re;
  result->c2.im += input->c2.im;
  result->c3.re += input->c3.re;
  result->c3.im += input->c3.im;

}

// result = s1 + s2
static inline void vector_add(su3_vector* result, const su3_vector* const s1,
                              const su3_vector* const s2)
{
  result->c1.re = s1->c1.re + s2->c1.re;
  result->c1.im = s1->c1.im + s2->c1.im;
  result->c2.re = s1->c2.re + s2->c2.re;
  result->c2.im = s1->c2.im + s2->c2.im;
  result->c3.re = s1->c3.re + s2->c3.re;
  result->c3.im = s1->c3.im + s2->c3.im;

}

// result = s1 + i* s2
static inline void vector_i_add(su3_vector* result, const su3_vector* const s1, 
                                const su3_vector* const s2)
{
  
  result->c1.re = s1->c1.re - s2->c1.im;
  result->c1.im = s1->c1.im + s2->c1.re;
  result->c2.re = s1->c2.re - s2->c2.im;
  result->c2.im = s1->c2.im + s2->c2.re;
  result->c3.re = s1->c3.re - s2->c3.im;
  result->c3.im = s1->c3.im + s2->c3.re;

}

// result -= input
static inline void vector_sub_assign(su3_vector* result, 
                                     const su3_vector* const input)
{

  result->c1.re -= input->c1.re;
  result->c1.im -= input->c1.im;
  result->c2.re -= input->c2.re;
  result->c2.im -= input->c2.im;
  result->c3.re -= input->c3.re;
  result->c3.im -= input->c3.im;

}

// result = s1 - i* s2
static inline void vector_i_sub(su3_vector* result, const su3_vector* const s1, 
                                const su3_vector* const s2)
{

  result->c1.re = s1->c1.re + s2->c1.im;
  result->c1.im = s1->c1.im - s2->c1.re;
  result->c2.re = s1->c2.re + s2->c2.im;
  result->c2.im = s1->c2.im - s2->c2.re;
  result->c3.re = s1->c3.re + s2->c3.im;
  result->c3.im = s1->c3.im - s2->c3.re;

}

// result = s1 - s2
static inline void vector_sub(su3_vector* result, const su3_vector* const s1,
                              const su3_vector* const s2)
{
  result->c1.re = s1->c1.re - s2->c1.re;
  result->c1.im = s1->c1.im - s2->c1.im;
  result->c2.re = s1->c2.re - s2->c2.re;
  result->c2.im = s1->c2.im - s2->c2.im;
  result->c3.re = s1->c3.re - s2->c3.re;
  result->c3.im = s1->c3.im - s2->c3.im;

}

// real part of scalar product of two su3 vector
// r = ((*input1) * (*input2)^*).re
static inline Float vector_prod_re(const su3_vector* const input1, 
                            const su3_vector* const input2)
{

  Float r;

  r = (input1->c1.re*input2->c1.re + input1->c1.im*input2->c1.im +
       input1->c2.re*input2->c2.re + input1->c2.im*input2->c2.im +
       input1->c3.re*input2->c3.re + input1->c3.im*input2->c3.im);

  return r;

}

// imaginary part of scalar product of two su3 vector
// r = ((*input1) * (*input2)^*).im
static inline Float vector_prod_im(const su3_vector* const input1,
                            const su3_vector* const input2)
{

  Float r;

  r = (input1->c1.re*input2->c1.im - input1->c1.im*input2->c1.re +
       input1->c2.re*input2->c2.im - input1->c2.im*input2->c2.re +
       input1->c3.re*input2->c3.im - input1->c3.im*input2->c3.re);

  return r;

}

// multiplication of SU(3) matrix and vector
// result = m * v
static inline void su3_multiply(su3_vector* result, const su3* const m,
                         const su3_vector* const v)
{

  result->c1.re = (m->c11.re*v->c1.re - m->c11.im*v->c1.im +
                   m->c12.re*v->c2.re - m->c12.im*v->c2.im +
                   m->c13.re*v->c3.re - m->c13.im*v->c3.im);
  result->c1.im = (m->c11.re*v->c1.im + m->c11.im*v->c1.re +
                   m->c12.re*v->c1.im + m->c12.im*v->c2.re +
                   m->c13.re*v->c1.im + m->c13.im*v->c3.re);

  result->c2.re = (m->c21.re*v->c1.re - m->c21.im*v->c1.im +
                   m->c22.re*v->c2.re - m->c22.im*v->c2.im +
                   m->c23.re*v->c3.re - m->c23.im*v->c3.im);
  result->c2.im = (m->c21.re*v->c1.im + m->c21.im*v->c1.re +
                   m->c22.re*v->c1.im + m->c22.im*v->c2.re +
                   m->c23.re*v->c1.im + m->c23.im*v->c3.re);

  result->c3.re = (m->c31.re*v->c1.re - m->c31.im*v->c1.im +
                   m->c32.re*v->c2.re - m->c32.im*v->c2.im +
                   m->c33.re*v->c3.re - m->c33.im*v->c3.im);
  result->c3.im = (m->c31.re*v->c1.im + m->c31.im*v->c1.re +
                   m->c32.re*v->c1.im + m->c32.im*v->c2.re +
                   m->c33.re*v->c1.im + m->c33.im*v->c3.re);

}

// multiplication of SU(3)^dagger and vector
// result = m^dagger * v
static inline void su3_dag_multiply(su3_vector* result, const su3* const m,
                             const su3_vector* const v)
{

  result->c1.re = (m->c11.re*v->c1.re + m->c11.im*v->c1.im +
                   m->c21.re*v->c2.re + m->c21.im*v->c2.im +
                   m->c31.re*v->c3.re + m->c31.im*v->c3.im);
  result->c1.im = (m->c11.re*v->c1.im - m->c11.im*v->c1.re +
                   m->c21.re*v->c2.im - m->c21.im*v->c2.re +
                   m->c31.re*v->c3.im - m->c31.im*v->c3.re);

  result->c2.re = (m->c12.re*v->c1.re + m->c12.im*v->c1.im +
                   m->c22.re*v->c2.re + m->c22.im*v->c2.im +
                   m->c32.re*v->c3.re + m->c32.im*v->c3.im);
  result->c2.im = (m->c12.re*v->c1.im - m->c12.im*v->c1.re +
                   m->c22.re*v->c2.im - m->c22.im*v->c2.re +
                   m->c32.re*v->c3.im - m->c32.im*v->c3.re);

  result->c3.re = (m->c13.re*v->c1.re + m->c13.im*v->c1.im +
                   m->c23.re*v->c2.re + m->c23.im*v->c2.im +
                   m->c33.re*v->c3.re + m->c33.im*v->c3.im);
  result->c3.im = (m->c13.re*v->c1.im - m->c13.im*v->c1.re +
                   m->c23.re*v->c2.im - m->c23.im*v->c2.re +
                   m->c33.re*v->c3.im - m->c33.im*v->c3.re);

}

//multiplication of SU(3) matrix and Dirac spinor
static inline void su3_spinor_mul(spinor* result, const su3* const m,
                           const spinor* const s)
{

  su3_multiply(&(result->c1), m, &(s->c1));
  su3_multiply(&(result->c2), m, &(s->c2));
  su3_multiply(&(result->c3), m, &(s->c3));
  su3_multiply(&(result->c4), m, &(s->c4));

}

// multiplication of SU(3)^dagger and Dirac spinor
static inline void su3_dag_spinor_mul(spinor* result, const su3* const m,
                               const spinor* const s)
{
  su3_dag_multiply(&(result->c1), m, &(s->c1));
  su3_dag_multiply(&(result->c2), m, &(s->c2));
  su3_dag_multiply(&(result->c3), m, &(s->c3));
  su3_dag_multiply(&(result->c4), m, &(s->c4));

}

// multiplication of real number and Dirac spinor
static inline void re_spinor_mul(spinor* result, const Float real,
                          const spinor* const in)
{

  re_vector_mul(&(result->c1), real, &(in->c1));
  re_vector_mul(&(result->c2), real, &(in->c2));
  re_vector_mul(&(result->c3), real, &(in->c3));
  re_vector_mul(&(result->c4), real, &(in->c4));

}

// multiplication of complex number and Dirac spinor
static inline void cm_spinor_mul(spinor* result, const complex z,
                                 const spinor* const in)
{

  cm_vector_mul(&(result->c1), z, &(in->c1));
  cm_vector_mul(&(result->c2), z, &(in->c2));
  cm_vector_mul(&(result->c3), z, &(in->c3));
  cm_vector_mul(&(result->c4), z, &(in->c4));

}

// += of two spinor
static inline void spinor_add_assign(spinor* result, spinor* in)
{

  vector_add_assign(&(result->c1), &(in->c1));
  vector_add_assign(&(result->c2), &(in->c2));
  vector_add_assign(&(result->c3), &(in->c3));
  vector_add_assign(&(result->c4), &(in->c4));

}

// addition of two spinor
static inline void spinor_add(spinor* result, spinor* s1, spinor* s2)
{
  vector_add(&(result->c1), &(s1->c1), &(s2->c2));
  vector_add(&(result->c2), &(s1->c2), &(s2->c2));
  vector_add(&(result->c3), &(s1->c3), &(s2->c3));
  vector_add(&(result->c4), &(s1->c4), &(s2->c4));

}

// -= of two spinor
static inline void spinor_sub_assign(spinor* result, spinor* in)
{

  vector_sub_assign(&(result->c1), &(in->c1));
  vector_sub_assign(&(result->c2), &(in->c2));
  vector_sub_assign(&(result->c3), &(in->c3));
  vector_sub_assign(&(result->c4), &(in->c4));
}

// subtraction of two spinor
static inline void spinor_sub(spinor* result, spinor* s1, spinor* s2)
{

  vector_sub(&(result->c1), &(s1->c1), &(s2->c2));
  vector_sub(&(result->c2), &(s1->c2), &(s2->c2));
  vector_sub(&(result->c3), &(s1->c3), &(s2->c3));
  vector_sub(&(result->c4), &(s1->c4), &(s2->c4));

}

// inner product of two spinor
static inline complex spinor_prod(spinor* s1, spinor* s2)
{

  complex z;
      
  Float x,y;
  x = (vector_prod_re(&(s1->c1), &(s2->c1)) +
       vector_prod_re(&(s1->c2), &(s2->c2)) +
       vector_prod_re(&(s1->c3), &(s2->c3)) +
       vector_prod_re(&(s1->c4), &(s2->c4)));

  y = (vector_prod_im(&(s1->c1), &(s2->c1)) +
       vector_prod_im(&(s1->c2), &(s2->c2)) +
       vector_prod_im(&(s1->c3), &(s2->c3)) +
       vector_prod_im(&(s1->c4), &(s2->c4)));
  
  z.re = x;
  z.im = y;

  return z;

}

// norm square of a spinor
static inline Float norm_square(spinor* s)
{

  Float real = 0.0;
  real =(vector_prod_re(&(s->c1), &(s->c1)) + 
         vector_prod_re(&(s->c2), &(s->c2))  + 
         vector_prod_re(&(s->c3), &(s->c3))  + 
         vector_prod_re(&(s->c4), &(s->c4)));

  return real;

}


#endif

