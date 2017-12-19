/**********************************************************************
 *
 * implimentation of routines in init.h
 *
 *********************************************************************/
#include<stdlib.h>
#include<stdio.h>
#include "init.h"

void init_gauge_unit(su3_field u)
{
  complex one;
  complex zero;
  one.re=1.0; one.im=0.0;
  zero.re=0.0; zero.im=0.0;

  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          for(int imu=0; imu<4; imu++)
            {

              (u[it][iz][iy][ix][imu]).c11 = one;
              (u[it][iz][iy][ix][imu]).c12 = zero;
              (u[it][iz][iy][ix][imu]).c13 = zero;
              (u[it][iz][iy][ix][imu]).c21 = zero;
              (u[it][iz][iy][ix][imu]).c22 = one;
              (u[it][iz][iy][ix][imu]).c23 = zero;
              (u[it][iz][iy][ix][imu]).c31 = zero;
              (u[it][iz][iy][ix][imu]).c32 = zero;
              (u[it][iz][iy][ix][imu]).c33 = one;

            }

}

void init_gauge_random(su3_field u, int iseed)
{
  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          for(int imu=0; imu<4; imu++)
            {
              int index = (it*NZ*NY*NX*4 + iz*NY*NX*4 
                           + iy*NX*4 + ix*4 + imu);
              rlx_init(0, iseed+index);

              random_su3(&(u[it][iz][iy][ix][imu]));

            }
}

void init_spinor_zero(spinor* s)
{
  complex c_zero = {0.0, 0.0};
  su3_vector zero = {c_zero, c_zero, c_zero};

  s->c1 = zero;
  s->c2 = zero;
  s->c3 = zero;
  s->c4 = zero;

}

void init_spinor_unit(spinor_field s)
{
  complex one;
  one.re=1.0; one.im=0.0;
  su3_vector unit;
  unit.c1 = unit.c2 = unit.c3 = one;

  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          {
            (s[it][iz][iy][ix]).c1 = unit;
            (s[it][iz][iy][ix]).c2 = unit;
            (s[it][iz][iy][ix]).c3 = unit;
            (s[it][iz][iy][ix]).c4 = unit;

          }

}

void init_spinor_random(spinor_field s, int iseed)
{

  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          {

            int index = (it*NZ*NY*NX + iz*NY*NX
                         + iy*NX + ix);
            rlx_init(0, iseed+index);

            random_spinor(&(s[it][iz][iy][ix]), 1.0);

          }
}

void su3_field_alloc(su3_field** u)
{
  
  (*u) = (su3_field*)calloc(NX*NY*NZ*NT*4*18, sizeof(Float));

}

void su3_field_free(su3_field** u)
{
  
  free(*u);
  *u = NULL;
}

void spinor_field_alloc(spinor_field** s)
{
 
  (*s) = (spinor_field*)calloc(NX*NY*NZ*NT*4*3*2, sizeof(Float));

}

void spinor_field_free(spinor_field** s)
{
 
  free(*s);
  *s = NULL;

}
