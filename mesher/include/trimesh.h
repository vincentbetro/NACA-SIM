#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "Linked_List.h"
#include "Point.h"
#include "Vector.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))

#ifndef trimesh_h
#define trimesh_h

// actual routine
int trimesh(const int npt, const int tdim, const int nb, int *nbs, int ***bs, Point pt[], int **tri);
// Wrapper function for use with individual coordinate arrays
int trimesh(const int npt, const int tdim, const int nb, int *nbs, int ***bs, double *x, double *y, int **tri);

#endif
