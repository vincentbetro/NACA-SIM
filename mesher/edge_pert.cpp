#include "edge_pert.h"

//takes node in question, list of naboring nodes (to make edges with), list of triangles and tri array for cost function, physical nodes and computational nodes (u,v) (as storage location for perturbations), current cost (to check for improvement) which as a pointer returns new cost...returns 0 or 1 for no change or change
int edge_pert(int nodein, List *triangles, List *edges, double **nodep, double **nodec, int **tri, double &current_cost)
  {
  int i, k; //counter
  Point p0, p1; //point placeholders
  Vector oppedge, uv; //opposite edge vectors, u and v vector
  double uvmag, cost; //magnitude of u and v, and cost at node 
  int test = 0; //init to no change

  //we loop thru all surrounding edges based on the node-node hash
  for (i = 0; i < edges->max; i++)
    {
    //use placeholder for less mem access
    k = edges->list[i];

    //we don't have edge with self
    if (k == nodein)
      continue;

    //now, determine edge unit vector
    p0 = Point(nodep[nodein][0],nodep[nodein][1]);
    p1 = Point(nodep[k][0],nodep[k][1]);
    oppedge = Vector(p0,p1);

    //normalize it to assure directionality preserved
    oppedge.normalize();

    //now, find the magnitude of the u and v by making it into a vector
    uv = Vector(nodec[nodein][0],nodec[nodein][1]);
    uvmag = uv.magnitude();

    //now, scale the unit normal with uvmag
    oppedge *= uvmag;

    //now, perturb node based on u and v (nodec) in normal direction
    nodep[nodein][0] += oppedge[0];
    nodep[nodein][1] += oppedge[1];

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
      nodep[nodein][0] -= oppedge[0];
      nodep[nodein][1] -= oppedge[1];
      }

    }

  //finally, return test
  return(test);
  }
