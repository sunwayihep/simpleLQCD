/******************************************************************
 *
 * initialization of gauge field and Dirac spinor field
 *
 *****************************************************************/
#ifndef __INIT_H_
#define __INIT_H_
#include "dirac.h"
#include "random.h"

extern void init_gauge_unit(su3_field u);

extern void init_gauge_random(su3_field u, int iseed);

extern void init_spinor_zero(spinor* s);

extern void init_spinor_unit(spinor_field s);

extern void init_spinor_random(spinor_field s, int iseed);

// allocate memory for su3_field
extern void su3_field_alloc(su3_field** u);

// deallocate memory for su3_field
extern void su3_field_free(su3_field** u);

// allocate memory for spinor_field
extern void spinor_field_alloc(spinor_field** s);

// deallocate memory for spinor_field
extern void spinor_field_free(spinor_field** s);

#endif
