#include "area_calc.h"

//takes empty cv area array and fills in, takes tag of boundaries, tri conn and physical nodes

void area_calc(double *cvareaG, int *tag, int **tri, int nt, double **nodep, int nn)
  {
  int n0, n1, n2, nd0, nd1, nd2; //node placeholders to avoid extra memory access
  int i, j, k, l; //counters
  Point corner1, corner2, corner3, corner4, p0, p1, p2, tricent; //point class placeholders, corners, centroid
  Vector u1, v1, norm, s0, s1, s2, norm0, norm1, norm2, norm0m, norm1m, norm2m; //vector class placeholders, side vectors, norms of side vectors, midpoint versions too
  Point mdpt0, mdpt1, mdpt2; //three edge midpoints for ease of use
 
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
          printf("\nYou have a triangle with more than three sides according to the loop counter in area_calc.cpp!\n");
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
       }
     }
   }

  return;
  } 
