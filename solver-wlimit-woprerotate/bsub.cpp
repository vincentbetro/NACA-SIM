#include "bsub.h"

void bsub(double **A, double *x, const int n)
  {
  int i, j; //counters

  //first, define U (DO NOT WRITE OVER A!)
  double **U;
  U = (double**)calloc(n,sizeof(double*));
  for (i = 0; i < n; i++)
    U[i] = (double*)calloc(n,sizeof(double));

  //now, init to zeros
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      U[i][j] = 0.0;         
  
  //now, fill in U with upper tri A
  for (i = 0; i < n; i++)
    {
    for (j = 0; j < n; j++)
      {
      if (i <= j)
        U[i][j] = A[i][j];
      }
    }

  //now, do backward subs
  for (i = n-1; i > -1; i--)
    {
    x[i] /= U[i][i];
    for (j = i-1; j > -1; j--)
      x[j] -= x[i]*U[j][i]; 
    }

  //clean up memory
  for (i = 0; i < n; i++)
    free(U[i]);
  freenull(U);

  return;
  } 
