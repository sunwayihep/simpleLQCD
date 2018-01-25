/***********************************************************************
 *
 * source file of Wilson-Dirac operator 
 * related routines
 *
 **********************************************************************/
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "dirac.h"
#include "init.h"


void one_pm_gamma_mul(spinor* adr, const int mu, const int pm,
                      const spinor* const in)
{

  if(pm == -1){
      switch(mu){
          // (1-gamma_0)* spinor
        case 0 :
            {
              vector_i_sub(&(adr->c1), &(in->c1), &(in->c4));
              vector_i_sub(&(adr->c2), &(in->c2), &(in->c3));
              I_vector_mul(&(adr->c3), &(adr->c2));
              I_vector_mul(&(adr->c4), &(adr->c1));

            };break;
          // (1-gamma_1)* spinor
        case 1 :
            {
              vector_add(&(adr->c1), &(in->c1), &(in->c4));
              vector_sub(&(adr->c2), &(in->c2), &(in->c3));
              re_vector_mul(&(adr->c3), -1.0, &(adr->c2));
              adr->c4 = adr->c1;

            };break;
          // (1-gamma_2)* spinor
        case 2 :
            {
              vector_i_sub(&(adr->c1), &(in->c1), &(in->c3));
              vector_i_add(&(adr->c2), &(in->c2), &(in->c4));
              I_vector_mul(&(adr->c3), &(adr->c1));
              mI_vector_mul(&(adr->c4), &(adr->c2));

            };break;
          // (1-gamma_3)* spinor
        case 3 :
            {
              vector_sub(&(adr->c1), &(in->c1), &(in->c3));
              vector_sub(&(adr->c2), &(in->c2), &(in->c4));
              re_vector_mul(&(adr->c3), -1.0, &(adr->c1));
              re_vector_mul(&(adr->c4), -1.0, &(adr->c2));

            };break;
        default:
            {
              printf("wrong parameter mu!");
              exit(1);
            }      
      }

  }else if(pm == 1){

      switch(mu){
          // (1+gamma_0)* spinor
        case 0 :
            {
              vector_i_add(&(adr->c1), &(in->c1), &(in->c4));
              vector_i_add(&(adr->c2), &(in->c2), &(in->c3));
              mI_vector_mul(&(adr->c3), &(adr->c2));
              mI_vector_mul(&(adr->c4), &(adr->c1));

            };break;
          // (1+gamma_1)* spinor
        case 1 :
            {
              vector_sub(&(adr->c1), &(in->c1), &(in->c4));
              vector_add(&(adr->c2), &(in->c2), &(in->c3));
              adr->c3 = adr->c2;
              re_vector_mul(&(adr->c4), -1.0, &(adr->c1));

            };break;
          // (1+gamma_2)* spinor
        case 2 :
            {
              vector_i_add(&(adr->c1), &(in->c1), &(in->c3));
              vector_i_sub(&(adr->c2), &(in->c2), &(in->c4));
              mI_vector_mul(&(adr->c3), &(adr->c1));
              I_vector_mul(&(adr->c4), &(adr->c2));

            };break;
          // (1+gamma_3)* spinor
        case 3 :
            {
              vector_add(&(adr->c1), &(in->c1), &(in->c3));
              vector_add(&(adr->c2), &(in->c2), &(in->c4));
              adr->c3 = adr->c1;
              adr->c4 = adr->c2;

            };break;
        default:
            {
              printf("wrong parameter mu!");
              exit(1);
            }      
      }

  }else{
      printf("wrong parameter pm!");
      exit(1);
  }

}

void dslash(spinor_field result, su3_field u, spinor_field input)
{
  spinor tmp1, tmp2, res;

  for(int it=0; it<NT; it++)
    for(int iz=0; iz<NZ; iz++)
      for(int iy=0; iy<NY; iy++)
        for(int ix=0; ix<NX; ix++)
          {

            init_spinor_zero(&res);

            // x-direction, mu=0
            // forward 
            su3_spinor_mul(&tmp1, &(u[it][iz][iy][ix][0]),
                           &(input[it][iz][iy][INC_X(ix)]));
            one_pm_gamma_mul(&tmp2, 0, -1, &tmp1);
            spinor_add_assign(&res, &tmp2);
            
            // backward 
            su3_dag_spinor_mul(&tmp1, &(u[it][iz][iy][DEC_X(ix)][0]),
                               &(input[it][iz][iy][DEC_X(ix)]));
            one_pm_gamma_mul(&tmp2, 0, 1, &tmp1);
            spinor_add_assign(&res, &tmp2);

            // y-direction, mu=1
            // forward
            su3_spinor_mul(&tmp1, &(u[it][iz][iy][ix][1]),
                           &(input[it][iz][INC_Y(iy)][ix]));
            one_pm_gamma_mul(&tmp2, 1, -1, &tmp1);
            spinor_add_assign(&res, &tmp2);
            
            // backward
            su3_dag_spinor_mul(&tmp1, &(u[it][iz][DEC_Y(iy)][ix][1]),
                               &(input[it][iz][DEC_Y(iy)][ix]));
            one_pm_gamma_mul(&tmp2, 1, 1, &tmp1);
            spinor_add_assign(&res, &tmp2);
            
            // z-direction, mu=2
            // forward
            su3_spinor_mul(&tmp1, &(u[it][iz][iy][ix][2]),
                           &(input[it][INC_Z(iz)][iy][ix]));
            one_pm_gamma_mul(&tmp2, 2, -1, &tmp1);
            spinor_add_assign(&res, &tmp2);
            
            //backward
            su3_dag_spinor_mul(&tmp1, &(u[it][DEC_Z(iz)][iy][ix][2]),
                               &(input[it][DEC_Z(iz)][iy][ix]));
            one_pm_gamma_mul(&tmp2, 2, 1, &tmp1);
            spinor_add_assign(&res, &tmp2);

            // t-direction, mu=3
            // forward
            su3_spinor_mul(&tmp1, &(u[it][iz][iy][ix][3]),
                           &(input[INC_T(it)][iz][iy][ix]));
            one_pm_gamma_mul(&tmp2, 3, -1, &tmp1);
            spinor_add_assign(&res, &tmp2);

            su3_dag_spinor_mul(&tmp1, &(u[DEC_T(it)][iz][iy][ix][3]),
                               &(input[DEC_T(it)][iz][iy][ix]));
            one_pm_gamma_mul(&tmp2, 3, 1, &tmp1);
            spinor_add_assign(&res, &tmp2);

            result[it][iz][iy][ix] = res;

          }
}

void m_wilson(spinor_field result, su3_field u, spinor_field input)
{
  Float c1 = 1/(2*kappa), c2 = 0.5;
  spinor_field tmp1, tmp2;

  re_spinor_field_mul(result, c1, input);

  dslash(tmp1, u, input);
  re_spinor_field_mul(tmp2, c2, tmp1);

  spinor_field_sub_assign(result, tmp2);
}


void mr_solver(spinor_field solution, su3_field u, spinor_field source)
{
  spinor_field  tmp, r, Mr;

  init_spinor_field_zero(solution);
  m_wilson(tmp, u, solution);
  spinor_field_sub(r, source, tmp);

  double residual = sqrt(norm_square_field(r));

  int count=0;

  complex coef1, coef2;
  double coef;

  while((count<max_iter) && (residual>tolerance))
    {
      ++count;

      m_wilson(Mr, u, r);

      coef1 = spinor_field_prod(Mr, r);
      coef = 1.0/(norm_square_field(Mr));
      coef *= omega;
      coef2 = re_complex_mul(coef1, coef);

      cm_spinor_field_mul(tmp, coef2, r); 
      spinor_field_add_assign(solution, tmp);

      cm_spinor_field_mul(tmp, coef2, Mr); 
      spinor_field_sub_assign(r, tmp);

      residual = sqrt(norm_square_field(r));

      printf("iteration: %d, residual: %.16e\n", count, residual);

    }

  printf("\n");
  printf("inversion converge after %d of iterations:\n",count);
  printf("last residual is %.16e\n\n", residual);

}





















