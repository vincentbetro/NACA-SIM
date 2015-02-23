#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
#define freenull(a) {if(a) free(a); a=NULL;}

#ifndef lin_solve_h
#define lin_solve_h

//lin_solve takes the number of internal loop iterations, number of nodes to loop thru, mat (rebuilt in mat_build), tag to avoid solving at boundries, u and v which are either physical nodes (winslow) or changes in location (L-E), convergence value to jump out at, relaxation factor, and compressed row storage index arrays
//it will pass out RMSerror for convergence check, and iterin by reference for L-E
double lin_solve(int &iterin, int nn, double ***mat, int *tag, double *u, double *v, double cvg, double omega, int *ia, int *ja, int *iau);

#endif 
