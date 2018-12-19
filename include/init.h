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

extern void init_spinor_field_zero(spinor_field s);

extern void init_spinor_field_unit(spinor_field s);

extern void init_spinor_field_random(spinor_field s, int iseed);

// allocate memory for su3_field
extern void su3_field_alloc(su3_field** u);

// deallocate memory for su3_field
extern void su3_field_free(su3_field** u);

// allocate memory for spinor_field
extern void spinor_field_alloc(spinor_field** s);

// deallocate memory for spinor_field
extern void spinor_field_free(spinor_field** s);

// read external gauge field
extern void read_gauge_field(su3_field* u, char filename[]);

// write fermion field to a file
// cause I don't have a propagator yet -_-
extern void write_ferm_field(spinor_field* s, char filename[], FILE* fp);

// swap endianess
extern void swap_endian(su3_field* u);

// set point source at space-time (t,z,y,x) and color=0, spin=0
// I don't want to deal with others whatever.
// Alright... I just deal with others...
extern void set_pt_source(spinor_field* s, int t, int z, int y, int x, int color, int spin);

#endif
