#include "Mesh.h"
#include "Lu4.h"
#include "fsub.h"
#include "bsub.h"

double Mesh::lin_solve(int &tempiter)
  {
  int i, j, k, flag; //counters and flag

  //set RMSerr high to assure we start through
  double RMSerr = 1000.0;

  //start iter at 0 each time through
  int iter = 0;

  //now, allocate for sol matrix (actually a vector 4*nn, but for ease of access, we make 2-D)
  //this is our delta q to be fed back in to change q at end of routine
  //resets to zero each time (as delta q is zero after each recompute)
  double **sol;
  sol = (double**)calloc(nn,sizeof(double*));
  for (i = 0; i < nn; i++)
    sol[i] = (double*)calloc(4,sizeof(double));

  //init temp matrix and sol vector for factorizer, forward/back sub routine purposes
  double **A;
  A = (double**)calloc(4,sizeof(double*));
  for (i = 0; i < 4; i++)
    A[i] = (double*)calloc(4,sizeof(double));

  double *x;
  x = (double*)calloc(4,sizeof(double));

  double Ax = 0.0; //will hold an element of vector Ax for finding RMS error cleanly

  //now, we perform our solve...we do this a prescribed number of times or until convergence
  do
    {
    iter++; //increment iteration counter

    //now, we loop through nodes and do solve...we are updating sol as we go, so this is gauss-seidel
    //in order to rectify bias situations, we will do the loop in the opposite direction every other time thru (symm GS-SSOR)  
    if (iter % 2 != 0)
      {
      for (i = 0; i < nn; i++)
        {
        //now, set x to be RHS from current node
        for (j = 0; j < 4; j++)
          x[j] = RHS[i][j];

        //now, transfer off-diag elements times sol to current node's "RHS" (x)
        for (j = ia[i]; j < ia[i+1]; j++)
          {
          if (ja[j] == i)
            continue; //this is the diagonal element
          //now, transfer each off diagonal times sol over changing signs and "summing"
          x[0] -= mat[j][0][0]*sol[ja[j]][0] + mat[j][0][1]*sol[ja[j]][1] + mat[j][0][2]*sol[ja[j]][2] + mat[j][0][3]*sol[ja[j]][3];
          x[1] -= mat[j][1][0]*sol[ja[j]][0] + mat[j][1][1]*sol[ja[j]][1] + mat[j][1][2]*sol[ja[j]][2] + mat[j][1][3]*sol[ja[j]][3];
          x[2] -= mat[j][2][0]*sol[ja[j]][0] + mat[j][2][1]*sol[ja[j]][1] + mat[j][2][2]*sol[ja[j]][2] + mat[j][2][3]*sol[ja[j]][3];
          x[3] -= mat[j][3][0]*sol[ja[j]][0] + mat[j][3][1]*sol[ja[j]][1] + mat[j][3][2]*sol[ja[j]][2] + mat[j][3][3]*sol[ja[j]][3];
          }

        //now, place mat[iau[i]] into temp matrix to send to factorizer (we need original diag matrices for RMS comp)
        for (j = 0; j < 4; j++)
          for (k = 0; k < 4; k++)
            A[j][k] = mat[iau[i]][j][k];

        //now, send 4x4 from main diag to LU factorizer
        Lu4(A,4);

        //now, do forward sub
        fsub(A, x, 4);

        //now, do backward sub
        bsub(A, x, 4);

        /*//direct method
        for (j = 0; j < 4; j++)
          x[j] /= mat[iau[i]][j][j];*/

        //now, update sol, relax delta
        for (j = 0; j < 4; j++)
          sol[i][j] += omega*(x[j] - sol[i][j]); 
        }
      }

    if (iter % 2 == 0)
      {
      for (i = nn-1; i > -1; i--)
        {
        //now, set x to be RHS from current node
        for (j = 0; j < 4; j++)
          x[j] = RHS[i][j];

        //now, transfer off-diag elements times sol to current node's "RHS" (x)
        for (j = ia[i]; j < ia[i+1]; j++)
          {
          if (ja[j] == i)
            continue; //this is the diagonal element
          //now, transfer each off diagonal times sol over changing signs and "summing"
          x[0] -= mat[j][0][0]*sol[ja[j]][0] + mat[j][0][1]*sol[ja[j]][1] + mat[j][0][2]*sol[ja[j]][2] + mat[j][0][3]*sol[ja[j]][3];
          x[1] -= mat[j][1][0]*sol[ja[j]][0] + mat[j][1][1]*sol[ja[j]][1] + mat[j][1][2]*sol[ja[j]][2] + mat[j][1][3]*sol[ja[j]][3];
          x[2] -= mat[j][2][0]*sol[ja[j]][0] + mat[j][2][1]*sol[ja[j]][1] + mat[j][2][2]*sol[ja[j]][2] + mat[j][2][3]*sol[ja[j]][3];
          x[3] -= mat[j][3][0]*sol[ja[j]][0] + mat[j][3][1]*sol[ja[j]][1] + mat[j][3][2]*sol[ja[j]][2] + mat[j][3][3]*sol[ja[j]][3];
          }

        //now, place mat[iau[i]] into temp matrix to send to factorizer (we need original diag matrices for RMS comp)
        for (j = 0; j < 4; j++)
          for (k = 0; k < 4; k++)
            A[j][k] = mat[iau[i]][j][k];

        //now, send 4x4 from main diag to LU factorizer
        Lu4(A,4);

        //now, do forward sub
        fsub(A, x, 4);

        //now, do backward sub
        bsub(A, x, 4);

        /*//direct method
        for (j = 0; j < 4; j++)
          x[j] /= mat[iau[i]][j][j];*/

        //now, update sol, relax delta
        for (j = 0; j < 4; j++)
          sol[i][j] += omega*(x[j] - sol[i][j]);  
        }
      }
     
    RMSerr = 0.0; //reset RMSerr to begin summing

    //now, find RMS error (b - Ax), using fully transferred RHS
    for (i = 0; i < nn; i++)
      {
      //now, set x to be orig RHS from current node
      for (j = 0; j < 4; j++)
        x[j] = RHS[i][j];

      //now, transfer off-diag elements times CURRENT sol to current node's "RHS" (x)
      for (j = ia[i]; j < ia[i+1]; j++)
        {
        if (ja[j] == i)
          continue; //this is the diagonal element
        //now, transfer each off diagonal times sol over changing signs and "summing"
        x[0] -= mat[j][0][0]*sol[ja[j]][0] + mat[j][0][1]*sol[ja[j]][1] + mat[j][0][2]*sol[ja[j]][2] + mat[j][0][3]*sol[ja[j]][3];
        x[1] -= mat[j][1][0]*sol[ja[j]][0] + mat[j][1][1]*sol[ja[j]][1] + mat[j][1][2]*sol[ja[j]][2] + mat[j][1][3]*sol[ja[j]][3];
        x[2] -= mat[j][2][0]*sol[ja[j]][0] + mat[j][2][1]*sol[ja[j]][1] + mat[j][2][2]*sol[ja[j]][2] + mat[j][2][3]*sol[ja[j]][3];
        x[3] -= mat[j][3][0]*sol[ja[j]][0] + mat[j][3][1]*sol[ja[j]][1] + mat[j][3][2]*sol[ja[j]][2] + mat[j][3][3]*sol[ja[j]][3];
        }

      //this will not zero initially, since RHS will not be everywhere 0        
      for (j = 0; j < 4; j++)
        {
        //first, discern Ax using current sol
        Ax = mat[iau[i]][j][0]*sol[i][0] + mat[iau[i]][j][1]*sol[i][1] + mat[iau[i]][j][2]*sol[i][2] + mat[iau[i]][j][3]*sol[i][3];
        //now compute RMS error for this row
        RMSerr += (x[j] - Ax)*(x[j] - Ax);
        }
      }

    //finally, compute RMSerr
    RMSerr = sqrt(RMSerr/(4*nn));

    //debug
    //printf("\nLin_solve.cpp: Current RMS error (b-Ax) of %16.10e after %d internal iterations\n",RMSerr,iter);
   
    }while (iter < tempiter && RMSerr > cvg);

  //set itermax as iter for passing out by ref
  tempiter = iter;

  //reset q
  for (i = 0; i < nn; i++)
    {
    for (j = 0; j < 4; j++)
      {
      q[i][j] += sol[i][j];
      }
    }

  //clean up memory
  for (i = 0; i < nn; i++)
    free(sol[i]);
  freenull(sol);
  freenull(A);
  freenull(x);

  //return RMS error
  return(RMSerr);
  } 
