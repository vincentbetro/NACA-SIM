#include "tri_edge.h"

//takes node in question, list of naboring triangles, physical nodes and computational nodes (u,v) (as storage location for perturbations), tri array, current cost (to check for improvement) which as a pointer returns new cost...returns 0 or 1 for no change or change
int tri_edge(int nodein, List *triangles, double **nodep, double **nodec, int **tri, double &current_cost)
  {
  int i, k; //counters
  int n0, n1, n2, nd0, nd1; //node placeholders
  Point p0, p1; //point placeholders
  Vector oppedge, oppnorm, uv; //opposite edge vectors, u and v vector
  double uvmag, cost; //magnitude of u and v, and cost at node
  int test = 0; //init to no change 

  //we loop thru all surrounding triangles based on the node-cell hash
  for (i = 0; i < triangles->max; i++)
    {
    //use placeholder for less mem access
    k = triangles->list[i]; 
    //list all vertices
    n0 = tri[k][0];
    n1 = tri[k][1];
    n2 = tri[k][2];
    //now, figure out which one is the node in question
    if (n0 == nodein)
      {
      nd0 = n1;
      nd1 = n2;
      }
    else if (n1 == nodein)
      {
      nd0 = n2;
      nd1 = n0;
      }
    else if (n2 == nodein)
      {
      nd0 = n0;
      nd1 = n1;
      }
    else
      {
      printf("\nYour triangle list doesn't include the node you are perturbing in tri_edge.cpp!\n");
      fflush(stdout);
      exit(0);
      }
    //now, determine opp edge normal
    p0 = Point(nodep[nd0][0],nodep[nd0][1]);
    p1 = Point(nodep[nd1][0],nodep[nd1][1]);
    oppedge = Vector(p0,p1);
    oppnorm = Vector(-oppedge[1],oppedge[0]);
    //normalize it to assure directionality preserved
    oppnorm.normalize();

    //now, find the magnitude of the u and v by making it into a vector
    uv = Vector(nodec[nodein][0],nodec[nodein][1]);
    uvmag = uv.magnitude();

    //now, scale the unit normal with uvmag
    oppnorm *= uvmag;

    //now, perturb node based on u and v (nodec) in normal direction
    nodep[nodein][0] += oppnorm[0];
    nodep[nodein][1] += oppnorm[1];

    //now, check and see if cost improved
    cost = cost_funct(triangles, tri, nodep);

    //now, determine if changes will stick
    if (cost < current_cost)
      {
      current_cost = cost; //reset cost
      //add perturbations to u and v
      //we increase by 10% but not ever decrease
      nodec[nodein][0] *= 1.1;
      nodec[nodein][1] *= 1.1;
      test = 1;
      }
    else
      {
      //unperturb node based on u and v (nodec) in normal direction
      nodep[nodein][0] -= oppnorm[0];
      nodep[nodein][1] -= oppnorm[1];
      }

    }

  //finally, return test
  return(test);
  }
