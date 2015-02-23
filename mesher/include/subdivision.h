#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Point.h"
#include "Vector.h"
#include "List.h"
#include "Linked_List.h"
#include "analytic_circ.h"
#include "analytic_curve.h"
#include "gradients.h"
#include "area_calc.h"
//eventually, we will include other functions that are functions of q

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
#define freenull(a) {if(a) free(a); a=NULL;}

//takes variable info, number of vars, number of functions to test on, which functions to test on, edge or grad based refinement toggle, tri connectivity, physical nodes, std deviation inputs for coarsening or refinement, node to node hash, nabor conn, boundry tagging, coarsening map, power to raise side length to for refinement function

//returns new node total

int subdivision(double **q, int nv, int test_funct, int **tri, int nt, double **nodep, int nn, int funct_type, double Cr[3], double Cc[3], double power[3], List **hash, int **nbr, int *tag, int *map, int ft); 
