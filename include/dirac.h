/**************************************************************************
 *
 * implementation of multiply of sparse Wilson-Dirac matrix and
 * spinor field 
 *
 * D_w = (N_d/a + m)*I - 1/2*Dslash
 * by defining kappa = 1/(2*(N_d/a+m))
 * D_w = 1/(2*kappa)*I - 1/2*Dslash
 *
 * where Dslash is
 * Dslash = \sum_\mu {(1-\gamma_\mu)U(x,\mu)\delta_{x+\mu,y} +
 *                   (1+\gamma_\mu)U(x-\mu,mu)^\dagger\delta_{x-\mu,y}}
 *
 * for Dirac gamma matrices, we choose chiral basis, the same as
 * QDP++ package, where gamma_5 is diagonal, 
 * gamma_0 = {{0,0,0,i},       
 *            {0,0,i,0},
 *            {0,-i,0,0},
 *            {-i,0,0,0}}
 *
 * gamma_1 = {{0,0,0,-1},
 *           {0,0,1,0},
 *           {0,1,0,0},
 *           {-1,0,0,0}}
 *
 * gamma_2 = {{0,0,i,0},
 *            {0,0,0,-i},
 *            {-i,0,0,0},
 *            {0,i,0,0}}
 *
 * gamma_3 = {{0,0,1,0},
 *            {0,0,0,1},
 *            {1,0,0,0},
 *            {0,1,0,0}}
 *
 * gamma_5 = {{1,0,0,0},
 *            {0,1,0,0},
 *            {0,0,-1,0},
 *            {0,0,0,-1}}
 *
 *************************************************************************/

#ifndef __DIRAC_H_
#define __DIRAC_H_

#include "linalg.h"

extern void one_pm_gamma_mul(spinor* add_to_result,int mu, int pm, 
                             const spinor* const input);

extern void dslash(spinor_field result, su3_field u, 
                   spinor_field input);

extern void m_wilson(spinor_field result, su3_field u, 
                     spinor_field input);

extern void mr_solver(spinor_field solution, su3_field u, 
                      spinor_field source);


#endif
