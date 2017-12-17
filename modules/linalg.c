/***************************************************
 *
 * source file of linear algebra routines
 *
 **************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "linalg.h"

complex spinor_field_prod(spinor_field s1, spinor_field s2)
{
  complex z;
  Float real=0.0, imag=0.0;

  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          {
            z = spinor_prod(&(s1[it][iz][iy][ix]),
                            &(s2[it][iz][iy][ix]));
            real += z.re;
            imag += z.im;
          
          }
  z.re = real;
  z.im = imag;

  return z;

}

Float norm_square_field(spinor_field s)
{
  Float r=0.0;

  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          {
            r += norm_square(&(s[it][iz][iy][ix]));
          
          }

  return r;
}

void re_spinor_field_mul(spinor_field result, const Float real,
                         spinor_field in)
{
  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          {
            re_spinor_mul(&(result[it][iz][iy][ix]), real,
                          &(in[it][iz][iy][ix]));
            
          }

}

void cm_spinor_field_mul(spinor_field result, const complex z,
                         spinor_field in)
{
  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          {
            cm_spinor_mul(&(result[it][iz][iy][ix]), z,
                          &(in[it][iz][iy][ix]));
            
          }

}

void spinor_field_add_assign(spinor_field result, spinor_field in)
{
  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          {
            spinor_add_assign(&(result[it][iz][iy][ix]),
                       &(in[it][iz][iy][ix]));
          
          }

}

void spinor_field_add(spinor_field result, spinor_field s1,
                      spinor_field s2)
{
  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          {
            spinor_add(&(result[it][iz][iy][ix]),
                       &(s1[it][iz][iy][ix]),
                       &(s2[it][iz][iy][ix]));
          
          }
}


void spinor_field_sub_assign(spinor_field result, spinor_field in)
{
  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          {
            spinor_sub_assign(&(result[it][iz][iy][ix]),
                       &(in[it][iz][iy][ix]));
          
          }

}

void spinor_field_sub(spinor_field result, spinor_field s1,
                      spinor_field s2)
{
  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          {
            spinor_sub(&(result[it][iz][iy][ix]),
                       &(s1[it][iz][iy][ix]),
                       &(s2[it][iz][iy][ix]));
          
          }
}
