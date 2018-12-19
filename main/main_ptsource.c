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

	su3_field* data;
	su3_field_alloc(&data);
	// test gauge configuration
	// 8x8x8x8 double precision
	char filename[] = "test.8.cfg";
	printf("size of su3_field = %i\n", sizeof(data));
	read_gauge_field(data, filename);
	// I need to swap endianess on my machine -_-
	swap_endian(data);

	printf("gauge field readed: \n");
	printf("%f\n", (((*data)[0][0][0][0][0]).c11.re));
	printf("%f\n", (((*data)[0][0][0][0][0]).c11.im));
	printf("%f\n", (((*data)[0][0][0][0][0]).c12.re));
	printf("%f\n", (((*data)[0][0][0][0][0]).c12.im));

	// write the solution from each spin and color source, i.e. 12 in total
	// to form a propagator
	char outPropName [] = "pt_source_prop_simpleLQCD";
	FILE* pf = fopen(outPropName, "w+b");

	spinor_field* check;
	spinor_field_alloc(&check);
	printf("Solving the Dirac equation to get propagator:\n");
	for(int spin=0; spin<4; ++spin)
		for(int color=0; color<3; ++color){
			set_pt_source(input, 0, 0, 0, 0, color, spin);
			printf("starting inversion: source spin= %d, source color= %d \n", spin, color);
			mr_solver(*solution, *data, *input);

			printf("Test for the correctness of MR inverter:\n\n");

			m_wilson(*check, *data, *solution);

			spinor_field_sub_assign(*check, *input);
			double res = sqrt(norm_square_field(*check));

			printf("True residual is:\n");
			printf("|M*solution - source| = %-20.16e\n\n",res);

		//	fwrite(solution, sizeof(spinor_field), 1, pf);
			write_ferm_field(solution, outPropName, pf);
		}
	fclose(pf);

	spinor_field_free(&check);
	su3_field_free(&data);
	spinor_field_free(&input);
	spinor_field_free(&solution);

	return 0;

}
