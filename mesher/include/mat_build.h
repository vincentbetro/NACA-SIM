#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Vector.h"
#include "Point.h"
#include "young_mod.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
#define freenull(a) {if(a) free(a); a=NULL;}

#ifndef mat_build_h
#define mat_build_h

//we pass in ux, uy, vx, vy to compute winslow coeff (these are sent in as 1 element arrays from L-E and ignored), poisson for L-E (set to zero for winslow), phi and psi (forcing functions), number of nodes, reinint to zeros version of mat to come back full, number of triangles, tri, from winslow, pass in nodec for node, from linear_elastic, pass in nodep for node, and compressed row storage indexing arrays
//the rub that won't let us compute these externally is that the alphas and betas are node based in winslow and the alphas and thetas are element based in L-E...this is why we must pass in smooth_type
//pass in weights for forcing functions and overall forcing function weight (C)
void mat_build(int smooth_type, double *ux, double *uy, double *vx, double *vy, double poisson, double *phi, double *psi, int nn, double ***mat, int nt, int **tri, double **node, int *ia, int *ja, int *iau, double *w1, double *w2, double *w1x, double *w1y, double *w2x, double *w2y, double C);

#endif 
