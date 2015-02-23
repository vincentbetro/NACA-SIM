#include "Mesh.h"

void Mesh::mesh_init(double rho0, double u0, double v0, double et0)
  {
  int n0, n1, n2, nd0, nd1, nd2; //node placeholders to avoid extra memory access
  int i, j, k, l; //counters
  Point corner1, corner2, corner3, corner4, p0, p1, p2, tricent; //point class placeholders, corners, centroid
  Vector u1, v1, norm, s0, s1, s2, norm0, norm1, norm2, norm0m, norm1m, norm2m; //vector class placeholders, side vectors, norms of side vectors, midpoint versions too
  Point mdpt0, mdpt1, mdpt2; //three edge midpoints for ease of use
  double mag0, mag1, mag2; //placeholders for magnitudes of vectors

  //to help later on, let's create tag array of boundary nodes
  tag = (int*)calloc(nn,sizeof(int));

  //init to ones
  for (i = 0; i < nn; i++)
    tag[i] = 1;

  //1 for interior, 0 for boundary
  for (i = 0; i < nb; i++)
    {
    for (j = 0; j < nbs[i]; j++)
      {
      for (k = 0; k < 2; k++)
        {
        if (tag[bs[i][j][k]] == 1)
          tag[bs[i][j][k]] = 0; //just to be sure we don't waste time overwriting
        }
      }
    }

  //to help later on, let's create tag array of boundary nodes, denoting inner (inviscid) or outer (free stream)
  tagbd = (int*)calloc(nn,sizeof(int));

  //init to zeroes
  for (i = 0; i < nn; i++)
    tagbd[i] = 0;

  //we only make some negative one if we have inner bd
  if (ib1 > -1 && ib2 > -1)
    {  
    //0 for outer boundary or non-boundary, 1 for inner boundary
    for (i = 0; i < 2; i++)
      {
      switch(i)
        {
        case 0:
          l = ib1;
        break;
        case 1:
          l = ib2;
        break;
        default:
          printf("\nYou chose a boundary that does not exist in tagbd switch in mesh_init.cpp!\n");
          fflush(stdout);
          exit(0);
        break;  
        }
      for (j = 0; j < nbs[l]; j++)
        {
        for (k = 0; k < 2; k++)
          {
          if (tagbd[bs[l][j][k]] == 0)
            tagbd[bs[l][j][k]] = 1; //just to be sure we don't waste time overwriting
          }
        }
      }
    }

  //allocate for cvareaG, len_scale
  cvareaG = (double*)calloc(nn,sizeof(double));
  len_scale = (double*)calloc(nn,sizeof(double));
  
  //we alloc another temp array to hold perimeter of each nodes cv
  double *perim = (double*)calloc(nn,sizeof(double));
    
  //now, we compute the cv areas using green's theorem
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

    //be careful about directions of normals for both of the following

    //now, we use green's theorem per triangle and sum areas associated with each node's cv, as well as take it away from the other node's cv as we go     
    cvareaG[n0] += 0.5*(norm0m*Vector(mdpt0[0],0.0) + norm0m*Vector(tricent[0],0.0));
    cvareaG[n1] -= 0.5*(norm0m*Vector(mdpt0[0],0.0) + norm0m*Vector(tricent[0],0.0));
    cvareaG[n0] -= 0.5*(norm2m*Vector(mdpt2[0],0.0) + norm2m*Vector(tricent[0],0.0));
    cvareaG[n2] += 0.5*(norm2m*Vector(mdpt2[0],0.0) + norm2m*Vector(tricent[0],0.0));
    cvareaG[n1] += 0.5*(norm1m*Vector(mdpt1[0],0.0) + norm1m*Vector(tricent[0],0.0));
    cvareaG[n2] -= 0.5*(norm1m*Vector(mdpt1[0],0.0) + norm1m*Vector(tricent[0],0.0));

    //begin to collectively sum each cv's perimeter
    mag0 = s0.magnitude();
    mag1 = s1.magnitude();
    mag2 = s2.magnitude();
    perim[n0] += mag0 + mag2;
    perim[n1] += mag0 + mag1;
    perim[n2] += mag1 + mag2; 

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
       //add in to cv area
       cvareaG[nd0] += 0.5*(norm0*Vector(corner1[0],0.0) + norm0*Vector(corner2[0],0.0));
       cvareaG[nd1] += 0.5*(norm1*Vector(corner3[0],0.0) + norm1*Vector(corner4[0],0.0));

       //add in to each cv's perimeter
       mag0 = s0.magnitude();
       mag1 = s1.magnitude();
       perim[nd0] += mag0;
       perim[nd1] += mag1;
       }
     }
   }


  //compute length scale by dividing each cv area by its perimeter
  for (i = 0; i < nn; i++)
    len_scale[i] = cvareaG[i]/perim[i];

  /*//debug
  for (i = 0; i < nn; i++)
    if (fabs(len_scale[i]) > 1.0e20 || fabs(len_scale[i]) < 1.0e-15)
       printf("\nlen_scale[%d] with tag = %d is bad.\n",i,tag[i]);*/

  //only set these if not read in
  if (nv == 0)
    {
    // allocate for each nodes q values
    q = (double**)calloc(nn,sizeof(double*));
    for (i = 0; i < nn; i++) 
      q[i] = (double*)calloc(4,sizeof(double));

    //set initial values of q
    for(i = 0; i < nn; i++)
      {
      q[i][0] = rho0;
      q[i][1] = rho0*u0;
      q[i][2] = rho0*v0;
      q[i][3] = et0;
      }
    }

  //clean up memory
  freenull(perim);

  return;
  } 
