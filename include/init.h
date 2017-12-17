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


#endif
