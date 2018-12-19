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

void init_spinor_field_zero(spinor_field s)
{

  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          {
          
            init_spinor_zero(&(s[it][iz][iy][ix]));
          
          }
}

void init_spinor_field_unit(spinor_field s)
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

void init_spinor_field_random(spinor_field s, int iseed)
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

void read_gauge_field(su3_field* u, char filename[]){
	FILE* ptr = fopen(filename, "rb");
	if(ptr == NULL){
		printf("Error in read configurations\n");
		exit(1);
	}
	fread(u, sizeof(su3_field), 1, ptr);
	fclose(ptr);
}

// tedious swap endianess function, but we're using C, so no complain-_-
void swap_endian(su3_field* u){
	const int size = sizeof(Float);
	char out[18][size];
	for(int it=0; it<NT; ++it)
		for(int iz=0; iz<NZ; ++iz)
			for(int iy=0; iy<NY; ++iy)
				for(int ix=0; ix<NX; ++ix)
					for(int imu=0; imu<4; ++imu)
					{
						for(int i=0; i<size; ++i){
							out[0][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c11.re)))[size-i-1];
							out[1][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c11.im)))[size-i-1];
							out[2][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c12.re)))[size-i-1];
							out[3][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c12.im)))[size-i-1];
							out[4][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c13.re)))[size-i-1];
							out[5][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c13.im)))[size-i-1];
							out[6][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c21.re)))[size-i-1];
							out[7][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c21.im)))[size-i-1];
							out[8][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c22.re)))[size-i-1];
							out[9][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c22.im)))[size-i-1];
							out[10][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c23.re)))[size-i-1];
							out[11][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c23.im)))[size-i-1];
							out[12][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c31.re)))[size-i-1];
							out[13][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c31.im)))[size-i-1];
							out[14][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c32.re)))[size-i-1];
							out[15][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c32.im)))[size-i-1];
							out[16][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c33.re)))[size-i-1];
							out[17][i] = ((char*)(&((*u)[it][iz][iy][ix][imu].c33.im)))[size-i-1];
						}


						for(int i=0; i<size; ++i){
							((char*)(&((*u)[it][iz][iy][ix][imu].c11.re)))[i] = out[0][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c11.im)))[i] = out[1][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c12.re)))[i] = out[2][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c12.im)))[i] = out[3][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c13.re)))[i] = out[4][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c13.im)))[i] = out[5][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c21.re)))[i] = out[6][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c21.im)))[i] = out[7][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c22.re)))[i] = out[8][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c22.im)))[i] = out[9][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c23.re)))[i] = out[10][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c23.im)))[i] = out[11][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c31.re)))[i] = out[12][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c31.im)))[i] = out[13][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c32.re)))[i] = out[14][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c32.im)))[i] = out[15][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c33.re)))[i] = out[16][i];
							((char*)(&((*u)[it][iz][iy][ix][imu].c33.im)))[i] = out[17][i];
						}
					}
}

void set_pt_source(spinor_field* s, int t, int z, int y, int x, int color, int spin){
	init_spinor_field_zero(*s);
	switch(spin){
		case 0:
			{
				switch(color){
					case 0:
						(*s)[t][z][y][x].c1.c1.re = 1.0;
						(*s)[t][z][y][x].c1.c1.im = 0.0;
						break;
					case 1:
						(*s)[t][z][y][x].c1.c2.re = 1.0;
						(*s)[t][z][y][x].c1.c2.im = 0.0;
						break;
					case 2:
						(*s)[t][z][y][x].c1.c3.re = 1.0;
						(*s)[t][z][y][x].c1.c3.im = 1.0;
						break;
				}
			}break;
		case 1:
			{
				switch(color){
					case 0:
						(*s)[t][z][y][x].c2.c1.re = 1.0;
						(*s)[t][z][y][x].c2.c1.im = 0.0;
						break;
					case 1:
						(*s)[t][z][y][x].c2.c2.re = 1.0;
						(*s)[t][z][y][x].c2.c2.im = 0.0;
						break;
					case 2:
						(*s)[t][z][y][x].c2.c3.re = 1.0;
						(*s)[t][z][y][x].c2.c3.im = 1.0;
						break;
				}
			}break;
		case 2:
			{
				switch(color){
					case 0:
						(*s)[t][z][y][x].c3.c1.re = 1.0;
						(*s)[t][z][y][x].c3.c1.im = 0.0;
						break;
					case 1:
						(*s)[t][z][y][x].c3.c2.re = 1.0;
						(*s)[t][z][y][x].c3.c2.im = 0.0;
						break;
					case 2:
						(*s)[t][z][y][x].c3.c3.re = 1.0;
						(*s)[t][z][y][x].c3.c3.im = 1.0;
						break;
				}
			}break;
		case 3:
			{
				switch(color){
					case 0:
						(*s)[t][z][y][x].c4.c1.re = 1.0;
						(*s)[t][z][y][x].c4.c1.im = 0.0;
						break;
					case 1:
						(*s)[t][z][y][x].c4.c2.re = 1.0;
						(*s)[t][z][y][x].c4.c2.im = 0.0;
						break;
					case 2:
						(*s)[t][z][y][x].c4.c3.re = 1.0;
						(*s)[t][z][y][x].c4.c3.im = 1.0;
						break;
				}
			}break;
		default:
			printf("Wrong source color and spin index, you need reconsider-_-\n");
			exit(1);
	}



}
