#include "opt_smooth.h"
#include "cost_funct.h"
#include "tri_edge.h"
#include "edge_pert.h"
#include "xy_pert.h"

//distance function to make finding epsilon easier...takes two points
double distance(Point p1, Point p2);

//takes tag, node-cell hash, node-node hash, physical nodes and computational nodes (u,v) (second only as storage location for perturbations), tri array, and dimensions into optimization-based smoothing procedure, iterations / cost at which to cut off convergence
void opt_smooth(int nn, double **nodep, double **nodec, List **hash, List **NChash, int nt, int **tri, int iterin, int *tag, double cvg)
  {
  int i, j, k; //counters
  Point p0, p1; //point placeholders
  double dst = 1.0e20; //distance between two pts, init high
  double eps = 1.0e20; //set epsilon high initially
  double cost = 0.0; //cost to make determination of smoothing options
  double min_cost, max_cost, avg_cost; //metrics that are tracked at end of each iteration
  int n_nodes; //number bd nodes for avg
  int test = 0; //output from functions to decide if x-y pert is necessary

  //first, we compute epsilon for the whole mesh by working through node hash and finding minimum (amongst all nodes) average distance between each node and its nabors...also, we create initial u and v here, then increment them as we go (and check for non-zero/below eps in the main loop)
  for (i = 0; i < nn; i++)
    {
    //compute u and v, but instead of u, v, we use nodec[0] and nodec[1]
    //we don't need u and v for bd, so they will be left 0.0
    nodec[i][0] = nodec[i][1] = 0.0; //init perturbations to zero
    //set up node placeholder for main node only once
    p0 = Point(nodep[i][0],nodep[i][1]);
    //loop thru naboring pts
    for (j = 0; j < hash[i]->max; j++)
      {
      //use placeholder to save mem access
      k = hash[i]->list[j];
      //we don't want distance between node and self
      if (k == i)
        continue; 
      //set up point placeholder for node's nabor
      p1 = Point(nodep[k][0],nodep[k][1]);
      //find distance
      dst = MAX(MIN(dst,distance(p0,p1)),1.0e-15);
      //set up u additively
      nodec[i][0] += nodep[k][0] - nodep[i][0];
      //set up v additively
      nodec[i][1] += nodep[k][1] - nodep[i][1];
      }
    //now, reset eps, using distance
    eps = MIN(eps,dst);
    //now, take avg, again noting we don't take distance from self
    nodec[i][0] /= (hash[i]->max - 1);
    nodec[i][1] /= (hash[i]->max - 1);
    }
  //finalize epsilon, be sure non-zero, mult by 1.0e-4 to make appropriately tolerant
  eps = MAX(1.0e-15,eps*1.0e-4);

  printf("\nEpsilon = %16.10e\n",eps);

  int iter = 0; //iteration counter

  //now, we perform smoothing iterations until iterin is reached
  do
    {
    iter++; //increment iteration counter

    //reset cost metrics
    max_cost = -1.0e20;
    min_cost = 1.0e20;
    avg_cost = 0.0;

    //reset num non-bd nodes
    n_nodes = 0;
    
    //now, we loop thru nodes, smoothing based on global eps and diff movement directions/techniques
    for (i = 0; i < nn; i++)
      {
      //if we have bd node, ignore it
      if (tag[i] == 0)
        continue;

      //increment number of non-bd nodes for avg
      n_nodes++;

      //reset test
      test = 0;

      //now, assure non-zero (take max of u/v and eps in abs val, then give sign of u/v)
      nodec[i][0] = SIGN(MAX(fabs(nodec[i][0]),eps),nodec[i][0]);
      nodec[i][1] = SIGN(MAX(fabs(nodec[i][1]),eps),nodec[i][1]);

      //compute current cost..pass in triangles containing node, tri, physical coords
      cost = cost_funct(NChash[i], tri, nodep);
     
      //notice: we actually move nodes as we go (instead of all at end)      
      //now, do smoothing with opposite edges (triangle)
      test += tri_edge(i, NChash[i], nodep, nodec, tri, cost);

      //now, do smoothing with edges 
      test += edge_pert(i, NChash[i], hash[i], nodep, nodec, tri, cost);

      //now, if neither has made a difference, we do x-y perturbations
      if (!test)
        test = xy_pert(i, NChash[i], nodep, nodec, tri, cost); 

      //now, track min, max, avg cost
      min_cost = MIN(min_cost,cost);
      max_cost = MAX(max_cost,cost);
      avg_cost += cost; //will divide at end
 
      } 
    //finally, divide avg_cost by number of nodes
    avg_cost /= n_nodes;

    //monitor cost functions and perturbations here
    printf("\nIteration %d, min cost = %16.10e, max cost = %16.10e, avg cost = %16.10e\n", iter, min_cost, max_cost, avg_cost);

    } while (iter < iterin && avg_cost > cvg);

  //now, we go back to main and check our mesh!
  return;
  }

//distance function to make finding epsilon easier...takes two points
double distance(Point p1, Point p2)
  {
  double dx = p2[0]-p1[0];
  double dy = p2[1]-p1[1];
  double mag = sqrt(dx*dx+dy*dy);
  return(mag);
  }
