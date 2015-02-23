#include "Mesh.h"

void Mesh::gradients()
  {
  int n0, n1, n2, nd0, nd1, nd2; //node placeholders to avoid extra memory access
  int i, j, k, l; //counters
  Point corner1, corner2, corner3, corner4, p0, p1, p2, tricent; //point class placeholders, corners, centroid
  Vector u1, v1, norm, s0, s1, s2, norm0, norm1, norm2, norm0m, norm1m, norm2m; //vector class placeholders, side vectors, norms of side vectors, midpoint versions too
  Point mdpt0, mdpt1, mdpt2; //three edge midpoints for ease of use

  //allocate for dqdx, dqdy (only need to check one, since if one, then both init)
  if (dqdx == 0)
    {
    dqdx = (double**)calloc(nn,sizeof(double*));
    dqdy = (double**)calloc(nn,sizeof(double*));
    for (i = 0; i < nn; i++)
      {
      dqdx[i] = (double*)calloc(4,sizeof(double));
      dqdy[i] = (double*)calloc(4,sizeof(double));
      }
    //printf("\nAllocating for dqdx, dqdy.");
    }
  else
    {   
    //simply re-init gradients to 0 (leave q alone)
    for (i = 0; i < nn; i++)
      {
      for (j = 0; j < 4; j++)
        {
        dqdx[i][j] = 0.0;
        dqdy[i][j] = 0.0;
        }
      }
    //printf("\nOnly redimensioning dqdx, dqdy.");
    }
    
  //limiter
  if (phi == 0)
    {
    phi = (double**)calloc(nn,sizeof(double*));
    for (i = 0; i < nn; i++)
      {
      phi[i] = (double*)calloc(4,sizeof(double));
      //init to 1.0
      for (j = 0; j < 4; j++)
        phi[i][j] = 1.0;
      }
    }
  else
    {   
    //simply re-init phi to 1.0
    for (i = 0; i < nn; i++)
      {
      for (j = 0; j < 4; j++)
        {
        phi[i][j] = 1.0;
        }
      }
    }

  //also, we need arrays to keep up with maxes and mins
if(qmax == 0)
  {
  qmax = (double**)calloc(nn,sizeof(double*));
  for (i = 0; i < nn; i++)
    {
    qmax[i] = (double*)calloc(4,sizeof(double));
    //init to qval
    for (j = 0; j < 4; j++)
      qmax[i][j] = q[i][j];
    }
  }
else
  {
  for (i = 0; i < nn; i++)
    {
    //init to qval
    for (j = 0; j < 4; j++)
      qmax[i][j] = q[i][j];
    }
  }

if (qmin == 0)
  {  
  qmin = (double**)calloc(nn,sizeof(double*));
  for (i = 0; i < nn; i++)
    {
    qmin[i] = (double*)calloc(4,sizeof(double));
    //init to qval
    for (j = 0; j < 4; j++)
      qmin[i][j] = q[i][j];
    }
  }
else
  {
  for (i = 0; i < nn; i++)
    {
    //init to qval
    for (j = 0; j < 4; j++)
      qmin[i][j] = q[i][j];
    }
  }

  //now, loop through edges to determine qmax and qmin (will get both directions for triangles, but if a boundary edge, need to do both now)
  for (i = 0; i < nt; i++)
    {
    //use node indices for cleanliness
    n0 = tri[i][0];
    n1 = tri[i][1];
    n2 = tri[i][2];
    
    for (j = 0; j < 4; j++)
      {
      qmax[n0][j] = MAX(qmax[n0][j],q[n1][j]);
      qmin[n0][j] = MIN(qmin[n0][j],q[n1][j]);
      qmax[n1][j] = MAX(qmax[n1][j],q[n2][j]);
      qmin[n1][j] = MIN(qmin[n1][j],q[n2][j]);
      qmax[n2][j] = MAX(qmax[n2][j],q[n0][j]);
      qmin[n2][j] = MIN(qmin[n2][j],q[n0][j]);
      
      //check for bd
      if (tag[n0] == 0 && tag[n1] == 0)
        {
        qmax[n1][j] = MAX(qmax[n1][j],q[n0][j]);
        qmin[n1][j] = MIN(qmin[n1][j],q[n0][j]);
        }
      if (tag[n1] == 0 && tag[n2] == 0)
        {
        qmax[n2][j] = MAX(qmax[n2][j],q[n1][j]);
        qmin[n2][j] = MIN(qmin[n2][j],q[n1][j]);
        }
      if (tag[n2] == 0 && tag[n0] == 0)
        {
        qmax[n0][j] = MAX(qmax[n0][j],q[n2][j]);
        qmin[n0][j] = MIN(qmin[n0][j],q[n2][j]);
        }
      }
    }

  //we allocate for temp arrays to hold extrapolated q values at mdpts and centroids
  double *qtricent = (double*)calloc(4,sizeof(double));
  double *qmdpt0 = (double*)calloc(4,sizeof(double));
  double *qmdpt1 = (double*)calloc(4,sizeof(double));
  double *qmdpt2 = (double*)calloc(4,sizeof(double));

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

    //determine q at midpoints and centroid
    for (k = 0; k < 4; k++)
      {
      qtricent[k] = (q[n0][k] + q[n1][k] + q[n2][k])/3.0;
      qmdpt0[k] = (q[n0][k] + q[n1][k])/2.0;
      qmdpt1[k] = (q[n1][k] + q[n2][k])/2.0;
      qmdpt2[k] = (q[n2][k] + q[n0][k])/2.0;
      }

    //now, we use green's theorem per triangle and gradients associated with each node, as well as take away from the other node's gradients as we go 
    for (k = 0; k < 4; k++)
      {    
      dqdx[n0][k] += 0.5*(norm0m*Vector(qmdpt0[k],0.0) + norm0m*Vector(qtricent[k],0.0));
      dqdy[n0][k] += 0.5*(norm0m*Vector(0.0,qmdpt0[k]) + norm0m*Vector(0.0,qtricent[k]));
      dqdx[n1][k] -= 0.5*(norm0m*Vector(qmdpt0[k],0.0) + norm0m*Vector(qtricent[k],0.0));
      dqdy[n1][k] -= 0.5*(norm0m*Vector(0.0,qmdpt0[k]) + norm0m*Vector(0.0,qtricent[k]));
    
      dqdx[n1][k] += 0.5*(norm1m*Vector(qmdpt1[k],0.0) + norm1m*Vector(qtricent[k],0.0));
      dqdy[n1][k] += 0.5*(norm1m*Vector(0.0,qmdpt1[k]) + norm1m*Vector(0.0,qtricent[k]));
      dqdx[n2][k] -= 0.5*(norm1m*Vector(qmdpt1[k],0.0) + norm1m*Vector(qtricent[k],0.0));
      dqdy[n2][k] -= 0.5*(norm1m*Vector(0.0,qmdpt1[k]) + norm1m*Vector(0.0,qtricent[k]));
    
      dqdx[n2][k] += 0.5*(norm2m*Vector(qmdpt2[k],0.0) + norm2m*Vector(qtricent[k],0.0));
      dqdy[n2][k] += 0.5*(norm2m*Vector(0.0,qmdpt2[k]) + norm2m*Vector(0.0,qtricent[k]));
      dqdx[n0][k] -= 0.5*(norm2m*Vector(qmdpt2[k],0.0) + norm2m*Vector(qtricent[k],0.0));
      dqdy[n0][k] -= 0.5*(norm2m*Vector(0.0,qmdpt2[k]) + norm2m*Vector(0.0,qtricent[k]));
      }    

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
          for (k = 0; k < 4; k++)
            {
            qmdpt0[k] = (q[n0][k] + q[n1][k])/2.0;
            }  
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
          for (k = 0; k < 4; k++)
            {
            qmdpt0[k] = (q[n1][k] + q[n2][k])/2.0;
            } 
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
          for (k = 0; k < 4; k++)
            {
            qmdpt0[k] = (q[n2][k] + q[n0][k])/2.0;
            } 
        break;
        default:
          printf("\nYou have a triangle with more than three sides according to the loop counter in mesh_init.cpp!\n");
          fflush(stdout);
          exit(0);
        break;
        }
       
     //check if we have a boundary node, do calculations on given side
     if (tag[nd0] == 0 && tag[nd1] == 0)
       {
       for (k = 0; k < 4; k++)
         {
         dqdx[nd0][k] += 0.5*(norm0*Vector(qmdpt0[k],0.0) + norm0*Vector(q[nd0][k],0.0));
         dqdy[nd0][k] += 0.5*(norm0*Vector(0.0,qmdpt0[k]) + norm0*Vector(0.0,q[nd0][k]));
         dqdx[nd1][k] += 0.5*(norm1*Vector(qmdpt0[k],0.0) + norm1*Vector(q[nd1][k],0.0));
         dqdy[nd1][k] += 0.5*(norm1*Vector(0.0,qmdpt0[k]) + norm1*Vector(0.0,q[nd1][k]));
         /*dqdx[nd0][k] += 0.75*(norm0*Vector(q[nd0][k],0.0)) + 0.25*(norm0*Vector(q[nd1][k],0.0));
         dqdy[nd0][k] += 0.75*(norm0*Vector(0.0,q[nd0][k])) + 0.25*(norm0*Vector(0.0,q[nd1][k]));
         dqdx[nd1][k] += 0.75*(norm1*Vector(q[nd1][k],0.0)) + 0.25*(norm1*Vector(q[nd0][k],0.0));
         dqdy[nd1][k] += 0.75*(norm1*Vector(0.0,q[nd1][k])) + 0.25*(norm1*Vector(0.0,q[nd0][k]));*/
         }
       }
     }
   }

  //divide all gradients by area of their CV
  for (i = 0; i < nn; i++)
    {
    for (j = 0; j < 4; j++)
      {
      dqdx[i][j] /= cvareaG[i];
      dqdy[i][j] /= cvareaG[i];
      }
    }


  //clean up memory
  freenull(qtricent);
  freenull(qmdpt0);
  freenull(qmdpt1);
  freenull(qmdpt2);

  return;
  } 
