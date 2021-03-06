#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
#define freenull(a) {if(a) free(a); a=NULL;}

void splitflux(double *ql, double *qr, double gamma, double xnorm, double ynorm, double rlen, double *fps, double *fms);
