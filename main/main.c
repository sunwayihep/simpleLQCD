/******************************************************
 *
 * main function to test Wilson-Dirac operator
 *
 *****************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "linalg.h"
#include "dirac.h"
#include "init.h"

double kappa = 0.005;
int max_iter = 1000;
double tolerance = 1.0e-8;
double omega = 0.2;

int main(int argc, char* argv[])
{
  su3_field u;
  spinor_field input;
  spinor_field solution;


  int nseed = 10;
  int istart = 1;

  for(int iseed=istart; iseed<nseed; iseed++)
    {

      printf("****************************************\n");

      printf("Randomly initialization of gauge field\n" 
             "and source vector for iseed=%d:\n\n",iseed);

      init_gauge_unit(u);
      init_spinor_unit(input);
      //init_gauge_random(u, iseed);
      //init_spinor_random(input, iseed);
      printf("Initialization finished!\n\n");

      printf("Starting solve the linear system:\n");
      mr_solver(solution, u, input);

      printf("Inversion finished!\n");
      printf("****************************************\n");
      // multiply back the solution with Wilson-Dirac operator
      // to test the MR solver

      printf("Test for the correctness of MR inverter:\n\n");
      spinor_field check;
      m_wilson(check, u, solution);

      spinor_field_sub_assign(check, input);
      double res = sqrt(norm_square_field(check));

      printf("Correctness check passed:\n");
      printf("|M*solution - source| = %-20.16e\n\n",res);


    }

  return 0;

}
