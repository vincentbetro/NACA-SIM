#include "winslow.h"

// takes number of nodes, node arrays, number of elements in mat, blank mat, ia, ja, iau, number of tri, tri, iteration limits, convergence value, omega, tag for boundaries
//poisson is passed in as 0.0 for consistency in mat_build 
//pass in A, B, C (forcing function coeff), whichff to tell which forcing function to create, and q to create forcing functions
void winslow(int nn, double **nodep, double **nodec, int mdim, double ***mat, int *ia, int *ja, int *iau, int nt, int **tri, int iterin, int iterout, double cvg, double omega, int *tag, double poisson, double A, double B, double C, int whichff, double **q)
  {
  int i, j, k; //counters

  double gamma = 1.4; //constant

  int smooth_type = 0; //this is based on which routine is called, we're in winslow, so it's 0

  //first, we must allocate for each node's coords, forcing functions, and actual function values, as well as areas around each point
  double *u = NULL, *v = NULL, *phi = NULL, *psi = NULL, *funct = NULL, *areapt = NULL, *garbage = NULL;

  u = (double*)calloc(nn,sizeof(double));
  v = (double*)calloc(nn,sizeof(double));
  phi = (double*)calloc(nn,sizeof(double)); //sent in as 0.0
  psi = (double*)calloc(nn,sizeof(double)); //sent in as 0.0
  funct = (double*)calloc(nn,sizeof(double));
  areapt = (double*)calloc(nn,sizeof(double));
  garbage = (double*)calloc(nn,sizeof(double));

  //u and v are simply our nodes' x and y coords
  for (i = 0; i < nn; i++)
    {
    u[i] = nodep[i][0];
    v[i] = nodep[i][1];
    }

  //now, we compute the forcing function based off of which was chosen (if C != 0.0)
  if (C != 0.0)
    {
    for (i = 0; i < nn; i++)
      {
      if (whichff == 0)
        funct[i] = (gamma - 1.0)*(q[i][3] - 0.5*((q[i][1]*q[i][1] + q[i][2]*q[i][2])/q[i][0]));
      else if (whichff == 1)
        funct[i] = sqrt((q[i][1]*q[i][1] + q[i][2]*q[i][2])/(q[i][0]*q[i][0]));
      else if (whichff == 2)
        {
        double pr = (gamma - 1.0)*(q[i][3] - 0.5*((q[i][1]*q[i][1] + q[i][2]*q[i][2])/q[i][0]));
        double ve = sqrt((q[i][1]*q[i][1] + q[i][2]*q[i][2])/(q[i][0]*q[i][0]));
        funct[i] = ve / sqrt((gamma*pr)/q[i][0]);
        }
      }

    //also, find areas around each point using area_calc
    area_calc(areapt, tag, tri, nt, nodep, nn);  
    }

  //we init ux, uy, vx, vy, fx, fy, fx2, fy2, w1, w2, w1x, w1y, w2x, w2y to null pointers
  double *ux = NULL, *vx = NULL, *uy = NULL, *vy = NULL, *fx = NULL, *fy = NULL, *fx2 = NULL, *fy2 = NULL, *w1 = NULL, *w2 = NULL, *w1x = NULL, *w1y = NULL, *w2x = NULL, *w2y = NULL;

  //if we are using winslow, we allocate for the gradients
  ux = (double*)calloc(nn,sizeof(double));
  uy = (double*)calloc(nn,sizeof(double));
  vx = (double*)calloc(nn,sizeof(double));
  vy = (double*)calloc(nn,sizeof(double));
  fx = (double*)calloc(nn,sizeof(double));
  fy = (double*)calloc(nn,sizeof(double));
  fx2 = (double*)calloc(nn,sizeof(double));
  fy2 = (double*)calloc(nn,sizeof(double));
  w1 = (double*)calloc(nn,sizeof(double));
  w2 = (double*)calloc(nn,sizeof(double));
  w1x = (double*)calloc(nn,sizeof(double));
  w1y = (double*)calloc(nn,sizeof(double));
  w2x = (double*)calloc(nn,sizeof(double));
  w2y = (double*)calloc(nn,sizeof(double));

  //now, since forcing functions don't change with smoothing, compute gradients here  (if C != 0)
  if (C != 0.0)
    {
    green_gauss(fx, fy, garbage, garbage, nn, nt, tri, nodec, funct, garbage);
    green_gauss(fx2, garbage, garbage, fy2, nn, nt, tri, nodec, fx, fy);

    //in order to normalize w1, w2, find max of fx, fy, fx2, fy2
    //init all low
    double maxfx = -1.0e20, maxfy = -1.0e20, maxfx2 = -1.0e20, maxfy2 = -1.0e20; 
    for (i = 0; i < nn; i++)
      {
      maxfx = MAX(ABS(fx[i]),ABS(maxfx));
      maxfy = MAX(ABS(fy[i]),ABS(maxfy));
      maxfx2 = MAX(ABS(fx2[i]),ABS(maxfx2));
      maxfy2 = MAX(ABS(fy2[i]),ABS(maxfy2));

      //maxfx = MAX(fx[i],maxfx);
      //maxfy = MAX(fy[i],maxfy);
      //maxfx2 = MAX(fx2[i],maxfx2);
      //maxfy2 = MAX(fy2[i],maxfy2);
      }

    //debug
    printf("\nMax function x gradient = %lf\n",maxfx);
    printf("\nMax function y gradient = %lf\n",maxfy);
    printf("\nMax gradient of function x gradient = %lf\n",maxfx2);
    printf("\nMax gradient of function y gradient = %lf\n",maxfy2);

    //now, compute w1, w2
    for (i = 0; i < nn; i++)
      {
      //w1[i] = 1.0 + A*areapt[i]*((fx[i]*fx[i])/(maxfx*maxfx)) + B*areapt[i]*((fx2[i]*fx2[i])/(maxfx2*maxfx2));
      //w2[i] = 1.0 + A*areapt[i]*((fy[i]*fy[i])/(maxfy*maxfy)) + B*areapt[i]*((fy2[i]*fy2[i])/(maxfy2*maxfy2));
      w1[i] = 1.0 + A*((fx[i]*fx[i])/(maxfx*maxfx)) + B*((fx2[i]*fx2[i])/(maxfx2*maxfx2));
      w2[i] = 1.0 + A*((fy[i]*fy[i])/(maxfy*maxfy)) + B*((fy2[i]*fy2[i])/(maxfy2*maxfy2));
      //w1[i] = 1.0 + A*areapt[i]*((fx[i]*fx[i])) + B*areapt[i]*((fx2[i]*fx2[i]));
      //w2[i] = 1.0 + A*areapt[i]*((fy[i]*fy[i])) + B*areapt[i]*((fy2[i]*fy2[i]));
      }

    //now, compute w1x, w1y, w2x, w2y 
    green_gauss(w1x, w1y, w2x, w2y, nn, nt, tri, nodec, w1, w2);
    }

  //now, we begin our outer loop
  //set iter to 0 to start iteration counter
  int iter = 0;
  //set RMSerr high to assure we start through
  double RMSerr = 1000.0;

  //set up temporary storage for iterin, since it will be passed out by reference and thus changed
  int itertemp = iterin;

  //keep looping until we reach convergence or exhaust the prescribed number of iterations
  do
    {
    iter++; //increment iteration counter

    //now, do green_gauss to compute coefficients..gradients passed in empty and out full
    green_gauss(ux, uy, vx, vy, nn, nt, tri, nodec, u, v);

    if (iter > 1)
      {
      //now, reset mat array (only do if not first iteration, since already done there by calloc)
      for (i = 0; i < mdim; i++)
        for (j = 0; j < 2; j++)
          for (k = 0; k < 2; k++)
            mat[i][j][k] = 0.0;
      }

    //now, rebuild mat
    mat_build(smooth_type, ux, uy, vx, vy, poisson, phi, psi, nn, mat, nt, tri, nodec, ia, ja, iau, w1, w2, w1x, w1y, w2x, w2y, C);

    //reset iterin each time since overwritten since passed by reference
    iterin = itertemp;

    //now, do linear solve, which returns current RMS error, and jumps out early if reached before iterin
    RMSerr = lin_solve(iterin, nn, mat, tag, u, v, cvg, omega, ia, ja, iau);
        
    }while (iter < iterout && RMSerr > cvg);

  printf("\nYou converged to RMS = %16.10e in %d external iterations of %d internal iterations each.\n",RMSerr,iter,itertemp);

  //finally, reset physical nodes with u and v to be passed back out, plotted, and sent to mesh file
  for (i = 0; i < nn; i++)
    {
    nodep[i][0] = u[i];
    nodep[i][1] = v[i];
    }
  
  //clean up memory
  freenull(u);
  freenull(v);
  freenull(ux);
  freenull(uy);
  freenull(vx);
  freenull(vy);
  freenull(fx);
  freenull(fy);
  freenull(fx2);
  freenull(fy2);
  freenull(w1);
  freenull(w2);
  freenull(w1x);
  freenull(w1y);
  freenull(w2x);
  freenull(w2y);
  freenull(phi);
  freenull(psi);
  freenull(funct);
  freenull(areapt);

  return;
  }
