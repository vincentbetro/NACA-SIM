#include "green_gauss.h"

//compute gradients of u and v (ux,uy,vx,vy) where passed in empty (to be reinit to 0)
//also, pass in nn for indexing, nt for looping over triangles, and tri, as well as nodec (since only used for winslow), and physical coords thru u and v
void green_gauss(double *ux, double *uy, double *vx, double *vy, int nn, int nt, int **tri, double **nodec, double *u, double *v)
  {
  int i; //counter

  double *area; //allocate for temp area array to hold nodal areas
  area = (double*)calloc(nn,sizeof(double)); 

  int n0, n1, n2; //node variables to keepm things clean and save on memory calls
  Vector u1, v1, norm0, norm1, norm2, s0, s1, s2; //vectors to compute triangle area, normals, sides (coresponding to opposite node)
  Point p0, p1, p2; //points to compose vectors 

  //re-initialize all passed in gradient arrays to zero
  for (i = 0; i < nn; i++)
    {
    ux[i] = 0.0;
    vx[i] = 0.0;
    uy[i] = 0.0;
    vy[i] = 0.0;
    }

  //now, go through and compute gradient across each triangle, and additively store this in node arrays, to divided by area summation (per node) later
  for (i = 0; i < nt; i++)
    {
    //use node indices for cleanliness
    n0 = tri[i][0];
    n1 = tri[i][1];
    n2 = tri[i][2];
      
    //create points and vectors for ease
    p0 = Point(nodec[n0][0],nodec[n0][1]);
    p1 = Point(nodec[n1][0],nodec[n1][1]);
    p2 = Point(nodec[n2][0],nodec[n2][1]);
    u1 = Vector(p1,p2);
    v1 = Vector(p1,p0);

    //store area in temp array for each node in question
    area[n0] += .5*(u1 % v1);
    area[n1] += .5*(u1 % v1);
    area[n2] += .5*(u1 % v1);
      
    //now, compute gradients to be additively stored 
    //create vectors of triangle sides (opposite node that is norm#)
    s0 = Vector(p1,p2);
    s1 = Vector(p2,p0);
    s2 = Vector(p0,p1);
    //create normal vectors
    norm0 = Vector(s0[1],-s0[0]);
    norm1 = Vector(s1[1],-s1[0]);
    norm2 = Vector(s2[1],-s2[0]);

    //now, additively store gradients across each element that touches node from triangle (areas cancel from dividing to get gradient across element and multiplying in wieghted average/green-gauss recovery of gradient for node)
    ux[n0] += (-u[n0]*norm0[0] - u[n1]*norm1[0] - u[n2]*norm2[0])/2.0;
    ux[n1] += (-u[n0]*norm0[0] - u[n1]*norm1[0] - u[n2]*norm2[0])/2.0;
    ux[n2] += (-u[n0]*norm0[0] - u[n1]*norm1[0] - u[n2]*norm2[0])/2.0;

    uy[n0] += (-u[n0]*norm0[1] - u[n1]*norm1[1] - u[n2]*norm2[1])/2.0;
    uy[n1] += (-u[n0]*norm0[1] - u[n1]*norm1[1] - u[n2]*norm2[1])/2.0;
    uy[n2] += (-u[n0]*norm0[1] - u[n1]*norm1[1] - u[n2]*norm2[1])/2.0;
  
    vx[n0] += (-v[n0]*norm0[0] - v[n1]*norm1[0] - v[n2]*norm2[0])/2.0;
    vx[n1] += (-v[n0]*norm0[0] - v[n1]*norm1[0] - v[n2]*norm2[0])/2.0;
    vx[n2] += (-v[n0]*norm0[0] - v[n1]*norm1[0] - v[n2]*norm2[0])/2.0;

    vy[n0] += (-v[n0]*norm0[1] - v[n1]*norm1[1] - v[n2]*norm2[1])/2.0;
    vy[n1] += (-v[n0]*norm0[1] - v[n1]*norm1[1] - v[n2]*norm2[1])/2.0;
    vy[n2] += (-v[n0]*norm0[1] - v[n1]*norm1[1] - v[n2]*norm2[1])/2.0;
    } 
    
  //now, go back and divide each gradient by area summation to have the final gradient at each node
  for (i = 0; i < nn; i++)
    {
    ux[i] /= area[i];
    uy[i] /= area[i];
    vx[i] /= area[i];
    vy[i] /= area[i];
    }

  /*//debug
    for (i = 0; i < nn; i++)
      {
      printf("\nux[%d] = %lf",i,ux[i]);
      printf("\nuy[%d] = %lf",i,uy[i]);
      printf("\nvx[%d] = %lf",i,vx[i]);
      printf("\nvy[%d] = %lf",i,vy[i]);
      }*/

  freenull(area); //clean up memory allocated in this routine

  return;
  } 
