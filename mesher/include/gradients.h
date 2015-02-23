#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Point.h"
#include "Vector.h" 

//computes grads of functions into empty dfdx, dfdy arrays, takes physical coords, tri conn, function values at each node

void gradients(double *cvareaG, int *tag, int **tri, int nt, double **nodep, int nn, double *funct, double *dfdx, double *dfdy); 
