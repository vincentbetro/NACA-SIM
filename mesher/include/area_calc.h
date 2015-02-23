#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Point.h"
#include "Vector.h"

//takes empty cv area array and fills in, takes tag of boundaries, tri conn and physical nodes 

void area_calc(double *cvareaG, int *tag, int **tri, int nt, double **nodep, int nn);
