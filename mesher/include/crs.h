#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "List.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
#define freenull(a) {if(a) free(a); a=NULL;}

#ifndef crs_h
#define crs_h

// hash takes number of nodes (dimension of List hash), number of elements in hash, hash table, and blank ia, ja, iau to be passed out full
void crs(int nn, int mdim, int *ia, int *ja, int *iau, List **hash);

#endif
