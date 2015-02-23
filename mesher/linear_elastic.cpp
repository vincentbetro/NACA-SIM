#include "linear_elastic.h"

// takes number of nodes, node arrays, number of elements in mat, blank mat, ia, ja, iau, number of tri, tri, iteration limits, convergence value, omega, tag for boundaries
// here, nodec is merely used as a place holder for boundary perturbations 
void linear_elastic(int nn, double **nodep, double **nodec, int mdim, double ***mat, int *ia, int *ja, int *iau, int nt, int **tri, int iterin, double cvg, double omega, int *tag, double poisson)
  {
  int i, j, k; //counters

  int smooth_type = 1; //this is based on which routine is called, we're in L-E, so it's 1

  //first, we must allocate for each node's coords and forcing functions
  double *u, *v, *phi, *psi;

  u = (double*)calloc(nn,sizeof(double));
  v = (double*)calloc(nn,sizeof(double));
  phi = (double*)calloc(nn,sizeof(double));
  psi = (double*)calloc(nn,sizeof(double));

  //since we have no forcing functions here, we init and leave phi and psi 0.0
  for (i = 0; i < nn; i++)
    phi[i] = psi[i] = 0.0;

  //u and v here are perturbations on the boundary, the interior should be initialized to 0.0
  for (i = 0; i < nn; i++)
    {
    if (tag[i] == 0)
      {
      u[i] = nodec[i][0] - nodep[i][0];
      v[i] = nodec[i][1] - nodep[i][1];
      }
    else
      {
      u[i] = 0.0; //to be safe, even though calloc was used
      v[i] = 0.0; //to be safe, even though calloc was used
      }
    }

  //we init ux, uy, vx, vy to null pointers
  double *ux = NULL, *vx = NULL, *uy = NULL, *vy = NULL;

  //since we must pass these unused to mat_build, we init to very small arrays
  ux = (double*)calloc(1,sizeof(double));
  uy = (double*)calloc(1,sizeof(double));
  vx = (double*)calloc(1,sizeof(double));
  vy = (double*)calloc(1,sizeof(double));

  //since we are not resetting mat based on new physical points (we iterate until RMS is reached with u and v as deltas, then we have our new mesh), we have no outer loop!
 
  //do reset mat each time, though, since we are moving the mesh
  for (i = 0; i < mdim; i++)
    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++)
        mat[i][j][k] = 0.0;

  //now, build mat (extra phi, psi are 0.0 placeholders since C = 0.0 by default here)
  mat_build(smooth_type, ux, uy, vx, vy, poisson, phi, psi, nn, mat, nt, tri, nodep, ia, ja, iau, phi, psi, phi, psi, phi, psi, 0.0);

  //now, do linear solve, which returns current RMS error, and jumps out early if reached before iterin
  //also, here we have it send back iter in place of iterin to tell us how many iterations it performed
  double RMSerr = lin_solve(iterin, nn, mat, tag, u, v, cvg, omega, ia, ja, iau);
    
  //finally, reset physical nodes with changes from u and v to affect the next round of mat_build
  //done for all nodes, including boundary since the original nodep on boundaries is NOT perturbed (and perturbation doesn't change since not solved for)
  for (i = 0; i < nn; i++)
    {  
    nodep[i][0] += u[i];
    nodep[i][1] += v[i];
    }

  //we get iterin by reference to be the value of iter from lin_solve
  printf("\nYou converged to RMS = %16.10e in %d iterations.\n",RMSerr, iterin);
 
  //clean up memory
  freenull(u);
  freenull(v);
  freenull(ux);
  freenull(uy);
  freenull(vx);
  freenull(vy);
  freenull(phi);
  freenull(psi);

  return;
  }
