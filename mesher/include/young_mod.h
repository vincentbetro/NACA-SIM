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

#ifndef young_mod_h
#define young_mod_h

//takes the corners of a triangle, computes aspect ratio and weighted condition number as options to use for Young's Modulus (stiffness keeps small cells from being overly squashed) 
double young_mod(Point p0, Point p1, Point p2);

#endif  
