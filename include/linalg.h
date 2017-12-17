/*************************************************************
 *
 * Dirac spinor field related linear algebra routines
 *
 * **********************************************************/

#ifndef __LINALG_H_
#define __LINALG_H_
#include "su3.h"

extern complex spinor_field_prod(spinor_field s1, spinor_field s2);
extern Float norm_square_field(spinor_field s);

extern void re_spinor_field_mul(spinor_field result, 
                                const Float real, spinor_field in);
extern void cm_spinor_field_mul(spinor_field result, 
                                const complex z, spinor_field in);

extern void spinor_field_add_assign(spinor_field result, spinor_field in);
extern void spinor_field_add(spinor_field result, spinor_field s1,
                             spinor_field s2);

extern void spinor_field_sub_assign(spinor_field result, spinor_field in);
extern void spinor_field_sub(spinor_field result, spinor_field s1,
                             spinor_field s2);
#endif
