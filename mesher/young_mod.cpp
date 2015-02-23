#include "young_mod.h"

int invert2X2(double m[2][2]); //use this function found below to make this easier!

//takes the corners of a triangle, computes aspect ratio and weighted condition number as options to use for Young's Modulus (stiffness keeps small cells from being overly squashed) 
double young_mod(Point p0, Point p1, Point p2)
  {
  //first, we will find aspect ratio, and divide it by area to give forgiveness to smaller cells
  
  //create vectors of sides
  Vector su1 = Vector(p1,p2); //also will be used in cond number, order ignored in mag used for aspect ratio
  Vector sv1 = Vector(p1,p0); //also will be used in cond number, order ignored in mag used for aspect ratio
  Vector sw1 = Vector(p0,p2);
  
  //find their magnitudes
  double len0 = su1.magnitude();
  double len1 = sv1.magnitude();
  double len2 = sw1.magnitude();

  //find the semiperimeter
  double semiper = 0.5*(len0+len1+len2);

  //find area
  double area_tri = 0.5*(su1 % sv1);

  //be sure not to divide by essentially zero area
  area_tri = MAX(area_tri,1.0e-15);

  //find aspect ratio
  double AR = (len0*len1*len2)/(8.0*(semiper-len0)*(semiper-len1)*(semiper-len2));

  //now, divide by area
  AR /= area_tri;

  //now, we find weighted condition number
  //we do this by finding the norm of A*inv(W) = [col([su1[0] su1[1]]) col([sv1[0] sv1[1]])] * [col([1 0]) col([-1/sqrt(3) 2/sqrt(3)])] and multiplying by norm of W*inv(A) = [col([1 0]) col([1/2 sqrt(3)/2])] * [col([sv1[1]/(su1[0]*sv1[1] - sv1[0]*su1[1]) -su1[1]/(su1[0]*sv1[1] - sv1[0]*su1[1])]) col([-sv1[0]/(su1[0]*sv1[1] - sv1[0]*su1[1]) su1[0]/(su1[0]*sv1[1] - sv1[0]*su1[1])])] then dividing by 2
  //of course, we use the froebinius norm, which in this case is simply the sum of the squares of the entries in each matrix mutliplication earlier
  
  //first, allocate for matrices (no need to delete since static allocation)
  double A[2][2], invA[2][2], W[2][2], invW[2][2];

  //now, store A and invert it
  invA[0][0]=A[0][0]=su1[0];
  invA[1][0]=A[1][0]=su1[1];
  invA[0][1]=A[0][1]=sv1[0];
  invA[1][1]=A[1][1]=sv1[1];

  //if not invertable, set CN very high, else it will place inv of A in invA
  if (!invert2X2(invA))
    return(1.0e20);

  //now, store W and invert it (we know it is invertable since static)
  invW[0][0]=W[0][0]=1.0;
  invW[1][0]=W[1][0]=0.0;
  invW[0][1]=W[0][1]=0.5;
  invW[1][1]=W[1][1]=sqrt(3.0)/2.0;
  invert2X2(invW);

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

  /*//matrix multiplication done by hand here, and broken down for ease of debugging
  //also, sum of squares of entries has been applied (norm taken)
  AinvW += su1[0]*su1[0];
  AinvW += (-su1[0]/sqrt(3.0) + 2.0*sv1[0]/sqrt(3.0))*(-su1[0]/sqrt(3.0) + 2.0*sv1[0]/sqrt(3.0));
  AinvW += su1[1]*su1[1];
  AinvW += (-su1[1]/sqrt(3.0) + 2.0*sv1[1]/sqrt(3.0))*(-su1[1]/sqrt(3.0) + 2.0*sv1[1]/sqrt(3.0));
  AinvW = sqrt(AinvW);
  
  WinvA += (sv1[1]/(su1[0]*sv1[1] - sv1[0]*su1[1])) - 0.5*(su1[1]/(su1[0]*sv1[1] - sv1[0]*su1[1]));
  WinvA += (-sv1[0]/(su1[0]*sv1[1] - sv1[0]*su1[1])) + 0.5*(su1[0]/(su1[0]*sv1[1] - sv1[0]*su1[1]));
  WinvA += -sqrt(3.0)*0.5*(su1[1]/(su1[0]*sv1[1] - sv1[0]*su1[1]));
  WinvA += sqrt(3.0)*0.5*(su1[0]/(su1[0]*sv1[1] - sv1[0]*su1[1]));
  WinvA = sqrt(WinvA);*/

  //now, find product of these two norms, and divide by 2 * area
  double CN = (AinvW * WinvA)/(2.0*area_tri);

  //now, we wish to return the max of these two metrics
  double ym = MAX(CN,AR);

  return(ym);
  }

//this routine will invert a 2x2 a la Steve Karman
int invert2X2(double m[2][2])
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
