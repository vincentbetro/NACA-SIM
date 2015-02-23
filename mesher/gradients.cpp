#include "gradients.h"

//computes grads of functions into empty dfdx, dfdy arrays, takes physical coords, tri conn, function values at each node

void gradients(double *cvareaG, int *tag, int **tri, int nt, double **nodep, int nn, double *funct, double *dfdx, double *dfdy)
  {
  int n0, n1, n2, nd0, nd1, nd2; //node placeholders to avoid extra memory access
  int i, j; //counters
  Point corner1, corner2, corner3, corner4, p0, p1, p2, tricent; //point class placeholders, corners, centroid
  Vector u1, v1, norm, s0, s1, s2, norm0, norm1, norm2, norm0m, norm1m, norm2m; //vector class placeholders, side vectors, norms of side vectors, midpoint versions too
  Point mdpt0, mdpt1, mdpt2; //three edge midpoints for ease of use

  //simply re-init gradients to 0
  for (i = 0; i < nn; i++)
    {
    dfdx[i] = 0.0;
    dfdy[i] = 0.0;
    }

  //we declare vars to hold extrapolated funct values at mdpts and centroids
  double ftricent = 0.0;
  double fmdpt0 = 0.0;
  double fmdpt1 = 0.0;
  double fmdpt2 = 0.0;

  //now, we compute the dqdx and dqdy using green's theorem
  for (i = 0; i < nt; i++)
    {
    //use node indices for cleanliness
    n0 = tri[i][0];
    n1 = tri[i][1];
    n2 = tri[i][2];
      
    //create point class placeholders for ease
    p0 = Point(nodep[n0][0],nodep[n0][1]);
    p1 = Point(nodep[n1][0],nodep[n1][1]);
    p2 = Point(nodep[n2][0],nodep[n2][1]);

    //create normal vectors from each side
    s0 = Vector(p0,p1);
    s1 = Vector(p1,p2);
    s2 = Vector(p2,p0);
    norm0 = Vector(s0[1],-s0[0]);
    norm1 = Vector(s1[1],-s1[0]);
    norm2 = Vector(s2[1],-s2[0]);

    //find triangle center by averaging the nodes
    tricent = (p0 + p1 + p2)/3.0;
    //printf("\n(%lf, %lf)\n",tricent[0],tricent[1]);

    //find midpoints
    mdpt0 = (p0 + p1)/2.0;
    mdpt1 = (p1 + p2)/2.0;
    mdpt2 = (p2 + p0)/2.0;

    //create normal vectors from each cv piece
    s0 = Vector(mdpt0,tricent);
    s1 = Vector(mdpt1,tricent);
    s2 = Vector(mdpt2,tricent);
    norm0m = Vector(s0[1],-s0[0]);
    norm1m = Vector(s1[1],-s1[0]);
    norm2m = Vector(s2[1],-s2[0]);

    //determine funct values at midpoints and centroid
    ftricent = (funct[n0] + funct[n1] + funct[n2])/3.0;
    fmdpt0 = (funct[n0] + funct[n1])/2.0;
    fmdpt1 = (funct[n1] + funct[n2])/2.0;
    fmdpt2 = (funct[n2] + funct[n0])/2.0;

    //now, we use green's theorem per triangle and gradients associated with each node, as well as take away from the other node's gradients as we go   
    dfdx[n0] += 0.5*(norm0m*Vector(fmdpt0,0.0) + norm0m*Vector(ftricent,0.0));
    dfdy[n0] += 0.5*(norm0m*Vector(0.0,fmdpt0) + norm0m*Vector(0.0,ftricent));
    dfdx[n1] -= 0.5*(norm0m*Vector(fmdpt0,0.0) + norm0m*Vector(ftricent,0.0));
    dfdy[n1] -= 0.5*(norm0m*Vector(0.0,fmdpt0) + norm0m*Vector(0.0,ftricent));
    
    dfdx[n1] += 0.5*(norm1m*Vector(fmdpt1,0.0) + norm1m*Vector(ftricent,0.0));
    dfdy[n1] += 0.5*(norm1m*Vector(0.0,fmdpt1) + norm1m*Vector(0.0,ftricent));
    dfdx[n2] -= 0.5*(norm1m*Vector(fmdpt1,0.0) + norm1m*Vector(ftricent,0.0));
    dfdy[n2] -= 0.5*(norm1m*Vector(0.0,fmdpt1) + norm1m*Vector(0.0,ftricent));
    
    dfdx[n2] += 0.5*(norm2m*Vector(fmdpt2,0.0) + norm2m*Vector(ftricent,0.0));
    dfdy[n2] += 0.5*(norm2m*Vector(0.0,fmdpt2) + norm2m*Vector(0.0,ftricent));
    dfdx[n0] -= 0.5*(norm2m*Vector(fmdpt2,0.0) + norm2m*Vector(ftricent,0.0));
    dfdy[n0] -= 0.5*(norm2m*Vector(0.0,fmdpt2) + norm2m*Vector(0.0,ftricent));    

    //now, we have to add in boundary contributions
    for (j = 0; j < 3; j++)
      {
      switch(j)
        {
        case 0:
          corner1 = p0;
          corner2 = mdpt0;
          corner3 = mdpt0;
          corner4 = p1;
          s0 = Vector(p0,mdpt0);
          s1 = Vector(mdpt0,p1);
          norm0 = Vector(s0[1],-s0[0]);
          norm1 = Vector(s1[1],-s1[0]);
          nd0 = n0;
          nd1 = n1;
          fmdpt0 = (funct[n0] + funct[n1])/2.0;  
        break;
        case 1:
          corner1 = p1;
          corner2 = mdpt1;
          corner3 = mdpt1;
          corner4 = p2;
          s0 = Vector(p1,mdpt1);
          s1 = Vector(mdpt1,p2);
          norm0 = Vector(s0[1],-s0[0]);
          norm1 = Vector(s1[1],-s1[0]);
          nd0 = n1;
          nd1 = n2;
          fmdpt0 = (funct[n1] + funct[n2])/2.0; 
        break;
        case 2:
          corner1 = p2;
          corner2 = mdpt2;
          corner3 = mdpt2;
          corner4 = p0;
          s0 = Vector(p2,mdpt2);
          s1 = Vector(mdpt2,p0);
          norm0 = Vector(s0[1],-s0[0]);
          norm1 = Vector(s1[1],-s1[0]);
          nd0 = n2;
          nd1 = n0;
          fmdpt0 = (funct[n2] + funct[n0])/2.0; 
        break;
        default:
          printf("\nYou have a triangle with more than three sides according to the loop counter in gradients.cpp!\n");
          fflush(stdout);
          exit(0);
        break;
        }
       
     //check if we have a boundary node, do calculations on given side
     if (tag[nd0] == 0 && tag[nd1] == 0)
       {
       dfdx[nd0] += 0.5*(norm0*Vector(fmdpt0,0.0) + norm0*Vector(funct[nd0],0.0));
       dfdy[nd0] += 0.5*(norm0*Vector(0.0,fmdpt0) + norm0*Vector(0.0,funct[nd0]));
       dfdx[nd1] += 0.5*(norm1*Vector(fmdpt0,0.0) + norm1*Vector(funct[nd1],0.0));
       dfdy[nd1] += 0.5*(norm1*Vector(0.0,fmdpt0) + norm1*Vector(0.0,funct[nd1]));
       /*dfdx[nd0] += 0.75*(norm0*Vector(funct[nd0],0.0)) + 0.25*(norm0*Vector(funct[nd1],0.0));
       dfdy[nd0] += 0.75*(norm0*Vector(0.0,funct[nd0])) + 0.25*(norm0*Vector(0.0,funct[nd1]));
       dfdx[nd1] += 0.75*(norm1*Vector(funct[nd1],0.0)) + 0.25*(norm1*Vector(funct[nd0],0.0));
       dfdy[nd1] += 0.75*(norm1*Vector(0.0,funct[nd1])) + 0.25*(norm1*Vector(0.0,funct[nd0]));*/
       }
     }
   }

  //divide all gradients by area of their CV
  for (i = 0; i < nn; i++)
    {
    dfdx[i] /= cvareaG[i];
    dfdy[i] /= cvareaG[i];
    }

  return;
  } 
