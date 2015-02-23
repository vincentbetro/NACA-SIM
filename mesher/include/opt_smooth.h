#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Vector.h"
#include "Point.h"
#include "List.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
#define freenull(a) {if(a) free(a); a=NULL;}

#ifndef opt_smooth_h
#define opt_smooth_h

//takes tag, node-cell hash, node-node hash, physical nodes and computational nodes (u,v) (second only as storage location for perturbations), tri array, and dimensions into optimization-based smoothing procedure, iterations / cost at which to cut off convergence
void opt_smooth(int nn, double **nodep, double **nodec, List **hash, List **NChash, int nt, int **tri, int iterin, int *tag, double cvg);

#endif 
