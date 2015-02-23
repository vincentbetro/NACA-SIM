#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Vector.h"
#include "Point.h"
#include "List.h"
#include "cost_funct.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
#define freenull(a) {if(a) free(a); a=NULL;}

#ifndef xy_pert_h
#define xy_pert_h

//takes node in question, list of naboring triangles and tri array for computing cost, physical nodes and computational nodes (u,v) (as storage location for perturbations), current cost (to check for improvement) which as a pointer returns new cost...returns 0 or 1 for no change or change
int xy_pert(int nodein, List *triangles, double **nodep, double **nodec, int **tri, double &current_cost);

#endif 
