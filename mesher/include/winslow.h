#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Vector.h"
#include "Point.h"
#include "green_gauss.h"
#include "mat_build.h"
#include "lin_solve.h"
#include "area_calc.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
#define freenull(a) {if(a) free(a); a=NULL;}

#ifndef winslow_h
#define winslow_h

// takes number of nodes, node arrays, number of elements in mat, blank mat, ia, ja, iau, number of tri, tri, iteration limits, convergence value, omega, tag for boundaries
//poisson is passed in as 0.0 for consistency in mat_build
//pass in A, B, C (forcing function coeff), whichff to tell which forcing function to create, and q to create forcing functions 
void winslow(int nn, double **nodep, double **nodec, int mdim, double ***mat, int *ia, int *ja, int *iau, int nt, int **tri, int iterin, int iterout, double cvg, double omega, int *tag, double poisson, double A, double B, double C, int whichff, double **q);

#endif
