#include "cost_funct.h"

int invert(double m[2][2]); //helps finding cond number
double cond_num(Point p0, Point p1, Point p2); //takes three corners of triangle and gives weighted condition number
double jac(Point p0, Point p1, Point p2); //takes three corners of triangle and gives avg jacobian

//takes list of triangles with that point in them, tri conn, physical nodes 
double cost_funct(List *triangles, int **tri, double **nodep)
  {
  int i, t; //counters
  int n0, n1, n2; //placeholders
  Point p0, p1, p2; //point placeholders

  double CN, JAC; //condition number and jacobian values

  double cmax = 0.0; //max cost over all triangles
  double cavg = 0.0; //avg cost (to be built upon, then averaged) over all triangles
  double cost = 0.0; //cost per triangle, or if blending, final cost
  double blend = 0.0; //blending factor

  //loop thru all triangles containing the point in question
  for (i = 0; i < triangles->max; i++)
    {
    //grab a triangle from the list and its indices
    t = triangles->list[i];
    n0 = tri[t][0];
    n1 = tri[t][1];
    n2 = tri[t][2];

    //grab physical points to pass to cond number and jacobian
    p0 = Point(nodep[n0][0],nodep[n0][1]);
    p1 = Point(nodep[n1][0],nodep[n1][1]);
    p2 = Point(nodep[n2][0],nodep[n2][1]);

    //find condition number (weighted)
    CN = cond_num(p0, p1, p2);
    
    //find average jacobian
    JAC = jac(p0, p1, p2);

    //now, determine cost for this triangle 
    if (JAC < 0)
      cost = 1.0 - JAC; //inverted element
    else
      cost = 1.0 - (1.0/CN);

    //now, square cost, in lieu of doing it in the actual above formulas
    cost *= cost;

    //start summing for avg
    cavg += cost;

    //check for cmax (always positive)
    cmax = MAX(cmax, cost);
    }
 
  //take average
  cavg /= triangles->max;

  //if our max cost is positive and high enough, return it
  if (cmax > 1.0)
    {
    return(cmax);
    }
  else
    {
    //now, we figure up blending factor
    blend = 0.5*(tanh((cmax - 1.0)*5.0) + 1.0);
    //now, we refigure new cost
    cost = blend*cmax + (1.0 - blend)*cavg;
    return(cost);
    }
  }

//this routine will invert a 2x2 a la Steve Karman
int invert(double m[2][2])
  {
  double a, b, c, d;
  a = m[0][0];
  b = m[0][1];
  c = m[1][0];
  d = m[1][1];
  double det = a*d-b*c;
  double tol = 1.0e-20;
  if (fabs(det) < tol)
    return(0);
  m[0][0] =  d/det;
  m[0][1] = -b/det;
  m[1][0] = -c/det;
  m[1][1] =  a/det;
  
  return(1);
  } 

//this routine gives cond number (weighted)
double cond_num(Point p0, Point p1, Point p2)
  {
  //create vectors of two sides
  Vector su1 = Vector(p1,p2); //also will be used in cond number, order ignored in mag used for aspect ratio
  Vector sv1 = Vector(p1,p0); //also will be used in cond number, order ignored in mag used for aspect ratio

  //we find weighted condition number by finding the norm of A*inv(W) = [col([su1[0] su1[1]]) col([sv1[0] sv1[1]])] * [col([1 0]) col([-1/sqrt(3) 2/sqrt(3)])] and multiplying by norm of W*inv(A) = [col([1 0]) col([1/2 sqrt(3)/2])] * [col([sv1[1]/(su1[0]*sv1[1] - sv1[0]*su1[1]) -su1[1]/(su1[0]*sv1[1] - sv1[0]*su1[1])]) col([-sv1[0]/(su1[0]*sv1[1] - sv1[0]*su1[1]) su1[0]/(su1[0]*sv1[1] - sv1[0]*su1[1])])] then dividing by 2
  //of course, we use the froebinius norm, which in this case is simply the sum of the squares of the entries in each matrix mutliplication earlier
  
  //first, allocate for matrices (no need to delete since static allocation)
  double A[2][2], invA[2][2], W[2][2], invW[2][2];

  //now, store A and invert it
  invA[0][0]=A[0][0]=su1[0];
  invA[1][0]=A[1][0]=su1[1];
  invA[0][1]=A[0][1]=sv1[0];
  invA[1][1]=A[1][1]=sv1[1];

  //if not invertable, set CN very high, else it will place inv of A in invA
  if (!invert(invA))
    return(1.0e20);

  //now, store W and invert it (we know it is invertable since static)
  invW[0][0]=W[0][0]=1.0;
  invW[1][0]=W[1][0]=0.0;
  invW[0][1]=W[0][1]=0.5;
  invW[1][1]=W[1][1]=sqrt(3.0)/2.0;
  invert(invW);

  //now, do matrix multiplication to get elements in two new matrices, square the results, and summatively store them (froebinius norm) in two dummy vars
  double AinvW = 0.0, WinvA = 0.0; //will eventually take sqrt at end

  AinvW += (A[0][0]*invW[0][0]+A[0][1]*invW[1][0])*(A[0][0]*invW[0][0]+A[0][1]*invW[1][0]);
  AinvW += (A[0][0]*invW[0][1]+A[0][1]*invW[1][1])*(A[0][0]*invW[0][1]+A[0][1]*invW[1][1]);
  AinvW += (A[1][0]*invW[0][0]+A[1][1]*invW[1][0])*(A[1][0]*invW[0][0]+A[1][1]*invW[1][0]);
  AinvW += (A[1][0]*invW[0][1]+A[1][1]*invW[1][1])*(A[1][0]*invW[0][1]+A[1][1]*invW[1][1]);
  AinvW = sqrt(AinvW);

  WinvA += (W[0][0]*invA[0][0]+W[0][1]*invA[1][0])*(W[0][0]*invA[0][0]+W[0][1]*invA[1][0]);
  WinvA += (W[0][0]*invA[0][1]+W[0][1]*invA[1][1])*(W[0][0]*invA[0][1]+W[0][1]*invA[1][1]);
  WinvA += (W[1][0]*invA[0][0]+W[1][1]*invA[1][0])*(W[1][0]*invA[0][0]+W[1][1]*invA[1][0]);
  WinvA += (W[1][0]*invA[0][1]+W[1][1]*invA[1][1])*(W[1][0]*invA[0][1]+W[1][1]*invA[1][1]);
  WinvA = sqrt(WinvA);

  //now, find product of these two norms, and divide by 2 (CANNOT USE AREA WEIGHTED, ELSE CAN GET CN < 1!!!!!)
  double CN = (AinvW * WinvA)/2.0;

  return(CN);
  }

//this routine gives avg jacobian
double jac(Point p0, Point p1, Point p2)
  {
  //create vectors of two sides
  Vector su1 = Vector(p1,p2); 
  Vector sv1 = Vector(p1,p0);

  //find jacobian
  double avgjac = su1 % sv1;

  return(avgjac);
  }
