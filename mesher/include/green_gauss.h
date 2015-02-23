#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Vector.h"
#include "Point.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
#define freenull(a) {if(a) free(a); a=NULL;}

#ifndef green_gauss_h
#define green_gauss_h

//compute gradients of u and v (ux,uy,vx,vy) where passed in empty (to be reinit to 0)
//also, pass in nn for indexing, nt for looping over triangles, and tri, as well as nodec (since only used for winslow), and physical coords thru u and v
void green_gauss(double *ux, double *uy, double *vx, double *vy, int nn, int nt, int **tri, double **nodec, double *u, double *v);

#endif 
