#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Vector.h"
#include "Point.h"
#include "mat_build.h"
#include "lin_solve.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
#define freenull(a) {if(a) free(a); a=NULL;}

#ifndef linear_elastic_h
#define linear_elastic_h

// takes number of nodes, node arrays, number of elements in mat, blank mat, ia, ja, iau, number of tri, tri, iteration limits, convergence value, omega, tag for boundaries 
// here, nodec is merely used as a place holder for boundary perturbations
void linear_elastic(int nn, double **nodep, double **nodec, int mdim, double ***mat, int *ia, int *ja, int *iau, int nt, int **tri, int iterin, double cvg, double omega, int *tag, double poisson);

#endif 
