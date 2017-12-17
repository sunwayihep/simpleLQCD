/********************************************************
 *
 * Luscher's random number generator
 *
 *******************************************************/
#ifndef __RANDOM_H_
#define __RANDOM_H_
#include "su3.h"

extern void ranlux(Float r[], int n);
extern void rlx_init(int level, int seed);

extern void gauss(Float r[], int n);

extern void random_spinor(spinor* s, Float sigma);
extern void random_su3_vector(su3_vector* v);
extern void random_su3(su3* u);

#endif
