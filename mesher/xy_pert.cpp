#include "xy_pert.h"

//takes node in question, list of naboring triangles and tri array for computing cost, physical nodes and computational nodes (u,v) (as storage location for perturbations), current cost (to check for improvement) which as a pointer returns new cost...returns 0 or 1 for no change or change
int xy_pert(int nodein, List *triangles, double **nodep, double **nodec, int **tri, double &current_cost)
  {
  double cost; //cost at node
  int test = 0; //init test to no change
  //we wish to do the same operation for u and v, and since we used nodec, we can run a loop (u first, v second)
  for (int i = 0; i < 2; i++)
    { 
    //first, perturb node by u
    nodep[nodein][i] += nodec[nodein][i];

    //now, check and see if cost improved
    cost = cost_funct(triangles, tri, nodep);

    //now, determine if changes will stick
    if (cost < current_cost)
      {
      current_cost = cost; //reset cost
      //add perturbations to u
      //we increase by 10%
      nodec[nodein][i] *= 1.1;
      test = 1; //set test as having changed
      }
    else
      {
      //unperturb node based on u (nodec) in normal direction
      nodep[nodein][i] -= nodec[nodein][i];
      //add perturbations to u
      //we decrease by 50%
      nodec[nodein][i] *= -0.5;
      }
    }

  //finally, return changed or not
  return(test);
  }
