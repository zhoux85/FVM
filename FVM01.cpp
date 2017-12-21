/*----------------------------------------------------------------------------
//    Example of FDM Method for 2D heat equation
//
//          u_t = u_{ xx } +u_{ yy } +f(x, t)
//
//    Test problem :
//      Exact solution : u(x, y, t) = exp(-t) cos(pi*x) cos(pi*y)
//      Source term : f(t, x, y) = -u + (2pi^2 - 1)exp(-t)cos(pi*x)cos(pi*y)
//      (x,y) is in (0,1) and (0,1)
-------------------------------------------------------------------------------*/

//////////////////////////////////////////////////////////////////////////
/* modified by Xiang Zhou, 2017/11/5 */

//#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>  
#include <sys/timeb.h>

#define PI 3.141592

int const nx = 10, ny = 10;//grid numbers
double dx = 1.0 / nx, dy = 1.0 / ny;//space step, 3cm*3cm
double D = 1;//D: diffusion coefficient cm^2/ms

/* Time Step */
double dt = 1e-05; // Time step (ms)
double t=0.0; // Time (ms)
int steps; // Number of Steps
int increment; // Loop Control Variable
double delta = 0.001;

double aw[nx + 1][ny + 1], ae[nx + 1][ny + 1], an[nx + 1][ny + 1], as[nx + 1][ny + 1], tem[nx + 1][ny + 1], tem_old[nx + 1][ny + 1],
ap_old[nx + 1][ny + 1], ap[nx + 1][ny + 1];

/* Voltage */
double U[nx + 2][ny + 2]; // Initial Voltage (mv)
double Unew[nx + 2][ny + 2]; // Initial Voltage (mv)
double U_old[nx + 2][ny + 2];
double Uerror[nx + 2][ny + 2]; //  error 

double A[nx + 1][ny + 1]; //  Ax=b 

double x[nx + 1], y[nx + 1];

FILE *single_ap;
void performance();
double uexact(double x, double y, double t);
double usource(double uu, double x, double y, double t);
void tam(double *A, double *d, double *x, int n);

int main(int argc, char* argv[])
{
	// (x,y)  coordinate position
	x[1] = 0.0;
	x[nx] = 1.0;
	int i;
	for (i = 2; i <= nx - 1; i++){
		x[i] = x[i - 1] + dx;
	}

	y[1] = 0.0;
	y[ny] = 1.0;
	int j;
	for (j = 2; j <= ny - 1; j++){
		y[j] = y[j - 1] + dy;
	}

	//compute the coefficients of discretization equations: aw, ae, as, an, ap_old, ap 
	for (i = 0; i <= nx + 1; i++){
		for (j = 0; j <= ny + 1; j++){
			aw[i][j] = D / dx;
			ae[i][j] = D / dx;
			as[i][j] = D / dy;
			an[i][j] = D / dy;
			ap_old[i][j] = dy*dx / dt;
		}
	}
	for (i = 1; i <= nx; i++){
		for (j = 1; j <= ny; j++){
			ap[i][j] = aw[i][j] + ae[i][j] + as[i][j] + an[i][j] + ap_old[i][j];
		}
	}

	/* Data File */
	FILE *fevaluation;
	fevaluation = fopen("fevaluation", "w");

	/* Time Loop Conditions */
	t = 0.0; // Time (ms)

	double U_old0[nx+2][ny+2];
	for (i = 0; i <= nx + 1; i++){
		for (j = 0; j <= ny + 1; j++){
			U[i][j] = uexact(i*dx, j*dy, 0); // Initial Voltage (mv)
			U_old[i][j] = U[i][j];
			U_old0[i][j] = U[i][j];
		}
	}

	struct timeb start, end;
	int diff;
	ftime(&start);

	double Uexact11[nx + 2][ny + 2] = { 0 };
	double dUdt[nx + 2][ny + 2] = { 0 };
	steps = 8 / dt;

	double belta[nx + 1];
	double b[nx + 1];
	double c[nx + 1];
	double a[nx + 1];

	double f[nx + 1][ny + 1];
	double y_temp[nx + 1];

	int ncount;
	double max_diff, tmp;
	for (ncount = 0; ncount <= steps; ncount++){
		for (i = 1; i <= nx; i++){
			//****no flux boundary conditions*****
			U[i][0] = U[i][1];
			U_old[i][0] = U[i][1];
			U[i][ny + 1] = U[i][ny];
			U_old[i][ny + 1] = U[i][ny];
		}
		for (j = 1; j <= ny; j++){
			U[0][j] = U[1][j];
			U_old[0][j] = U[0][j];
			U[nx + 1][j] = U[nx][j];
			U_old[nx + 1][j] = U[nx + 1][j];
		}

		//*********** step 1 *******
		do{
			for (j = 1; j <= ny; j++){
				//set column f: Ax=f
				for (i = 1; i <= nx; i++){
					if (i == 1){
						f[i][j] = aw[i][j] * U[i][j - 1] + as[i][j] * U[i - 1][j]
							+ ae[i][j] * U_old[i][j + 1] + ap[i][j] * U_old0[i][j];
					}
					else if (i == nx){
						f[i][j] = aw[i][j] * U[i][j - 1] + an[i][j] * U[i + 1][j]
							+ ae[i][j] * U_old[i][j + 1] + ap[i][j] * U_old0[i][j];
					}
					else{
						f[i][j] = aw[i][j] * U[i][j - 1] + ae[i][j] * U_old[i][j + 1]
							+ ap[i][j] * U_old0[i][j];
					}
				}

				//set matrix A: Ax=f
				for (i = 1; i <= nx; i++){
					b[i] = ap[i][j];
				}
				c[nx] = 0.0;
				a[1] = 0.0;
				for (i = 1; i <= nx - 1; i++){
					c[i] = -an[i][j];
					a[i + 1] = -as[i+1][j];
				}

				//TDMA: solve Ax=f
				belta[1] = c[1] / b[1];
				y_temp[1] = f[1][j] / b[1];
				for (i = 2; i < nx; i++){ //i = 2,3,...,n-1
					belta[i] = c[i] / (b[i] - a[i] * belta[i - 1]);
					y_temp[i] = (f[i][j] - a[i] * y_temp[i - 1]) / (b[i] - a[i] * belta[i - 1]);
				}
				y_temp[nx] = (f[nx][j] - a[nx] * y_temp[nx - 1]) / (b[nx] - a[nx] * belta[nx - 1]);
				U[nx][j] = y_temp[nx];
				for (i = nx - 1; i >= 1; i--){
					U[i][j] = y_temp[i] - belta[i] * U[i + 1][j];
				}
			}

			max_diff = 0.0;
			for (i = 1; i <= nx; i++){
				for (j = 1; j <= ny; j++){
					tmp = fabs(U[i][j] - U_old[i][j]);
					if (tmp>max_diff){
						max_diff = tmp;
					}
				}
			}
			//****no flux boundary conditions*****
			for (i = 1; i <= nx; i++){
				U[i][0] = U[i][1];
				U[i][ny + 1] = U[i][ny];
			}
			for (j = 1; j <= ny; j++){
				U[0][j] = U[1][j];
				U[nx + 1][j] = U[nx][j];
			}
			for (i = 0; i <= nx + 1; i++){
				for (j = 0; j <= ny + 1; j++){
					U_old[i][j] = U[i][j];
				}
			}
		} while (max_diff>1e-6);


		//*********** step 2 *******
		for (i = 1; i <= nx; i++){
			for (j = 1; j <= ny; j++){
				U[i][j] = U[i][j] + dt*usource(U[i][j], i*dx, j*dy, ncount*dt);
			}
		}
		for (i = 1; i <= nx; i++){
			for (j = 1; j <= ny; j++){
				U_old0[i][j] = U[i][j];
				U_old[i][j] = U[i][j];
			}
		}
		//*********** step 2 *******
		for (i = 0; i <= nx + 1; i++){
			for (j = 0; j <= ny + 1; j++){
				Uexact11[i][j] = uexact(i*dx, j*dy, ncount*dt); // Initial Voltage (mv)
			}
		}
		//for (i = 1; i <= nx; i++){
		//	for (j = 1; j <= ny; j++){
		//		printf("%g\t%g\t%g\n", i, j, U[i][j]);
		//	}
		//}
	}

	//four corners of the square
	U[0][0] = U[0][1];
	U[0][ny + 1] = U[0][ny];
	U[nx + 1][0] = U[nx][1];
	U[nx + 1][ny + 1] = U[nx + 1][ny];

	for (i = 0; i <= nx + 1; i++){
		for (j = 0; j <= ny + 1; j++){
			Uexact11[i][j] = uexact(i*dx, j*dy, steps*dt); // Initial Voltage (mv)
		}
	}

	for (i = 0; i <= nx + 1; i++){
		for (j = 0; j <= ny + 1; j++){
			Uerror[i][j] = fabs(Uexact11[i][j] - U[i][j]); // Initial Voltage (mv)
		}
	}

	//matrix norm L_1 error
	double sum_err = 0, err_L1 = 0;
	for (j = 0; j <= ny + 1; j++){
		for (i = 0; i <= nx + 1; i++){
			sum_err += Uerror[i][j];
		}
	}
	err_L1 = sum_err*dx*dy;

	ftime(&end);
	diff = (int)(1000.0*(end.time - start.time) + (end.millitm - start.millitm));
	fprintf(fevaluation, "%d s\n%g\n", diff / 1000, err_L1);
	fclose(fevaluation);
}

double uexact(double x, double y, double t){
	return  exp(-t)*cos(PI*x)*cos(PI*y);
}

double usource(double uu, double x, double y, double t){
	return -uu + (2 * PI*PI)*exp(-t)*cos(PI*x)*cos(PI*y);
}


