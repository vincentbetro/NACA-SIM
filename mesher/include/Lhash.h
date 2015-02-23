#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "List.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
#define freenull(a) {if(a) free(a); a=NULL;}

#ifndef Lhash_h
#define Lhash_h

// hash takes number of nodes (dimension of List Lhash), number of triangles, tri array to use to create hash table, and blank hash table to be passed out full
void Lhash(int nn, int nt, int **tri, List **hash);

#endif
