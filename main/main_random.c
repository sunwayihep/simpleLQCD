/******************************************************
 *
 * main function to test Wilson-Dirac operator
 *
 *****************************************************/
#include<byteswap.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "linalg.h"
#include "dirac.h"
#include "init.h"

double kappa = 0.135;
int max_iter = 1000;
double tolerance = 1.0e-12;
double omega = 1.1;

// check the endianess of the machine
int isBigEndian(){
	int test = 1;
	char* p = (char*)&test;
	return p[0] == 0;
}

int main(int argc, char* argv[])
{

  if(isBigEndian()){
	  printf("This system is Big Endian.\n");
  }else{
	  printf("This system is Little Endian.\n");
  }

  printf("kappa= %f, omega= %f\n",kappa, omega);


  spinor_field* input;
  spinor_field* solution;
  spinor_field_alloc(&input);
  spinor_field_alloc(&solution);

  su3_field* u;
  su3_field_alloc(&u);

  int nseed = 10;
  int istart = 1;

  for(int iseed=istart; iseed<nseed; iseed++)
    {

      printf("****************************************\n");

      printf("Randomly initialization of gauge field\n" 
             "and source vector for iseed=%d:\n\n",iseed);

//      init_gauge_unit(*u);
//      init_spinor_field_unit(*input);
      init_gauge_random(*u, iseed);
      init_spinor_field_random(*input, iseed);


      printf("Initialization finished!\n\n");

      printf("Starting solve the linear system:\n");
      mr_solver(*solution, *u, *input);

      printf("Inversion finished!\n");
      printf("****************************************\n");
      // multiply back the solution with Wilson-Dirac operator
      // to test the MR solver

      printf("Test for the correctness of MR inverter:\n\n");
      //  spinor_field check;
      spinor_field* check;
      spinor_field_alloc(&check);

      m_wilson(*check, *u, *solution);

      spinor_field_sub_assign(*check, *input);
      double res = sqrt(norm_square_field(*check));

      printf("True residual is:\n");
      printf("|M*solution - source| = %-20.16e\n\n",res);

      spinor_field_free(&check);

    }

  su3_field_free(&u);
  spinor_field_free(&input);
  spinor_field_free(&solution);

  return 0;

}
