/********************************************************
 *
 * global parameters for lattice simulation
 *
 *******************************************************/
#ifndef __GLOBAL_H_
#define __GLOBAL_H_

//base datatype
#define Float double

#define NX 12
#define NY 12
#define NZ 12
#define NT 144

// parameters related to bare quark mass 
// of Wilson-Dirac operator
extern double kappa;

// parameters for Minimal-Residual inverter
extern int max_iter;
extern double tolerance;
extern double omega;


#define INC_X(a) (((a)+1) % NX)
#define INC_Y(a) (((a)+1) % NY)
#define INC_Z(a) (((a)+1) % NZ)
#define INC_T(a) (((a)+1) % NT)
#define DEC_X(a) (((a)+NX-1) % NX)
#define DEC_Y(a) (((a)+NY-1) % NY)
#define DEC_Z(a) (((a)+NZ-1) % NZ)
#define DEC_T(a) (((a)+NT-1) % NT)


#endif
