#include "fsub.h"

void fsub(double **A, double *x, const int n)
  {

  int i, j; //counters

  //first, define L (DO NOT WRITE OVER A!)
  double **L;
  L = (double**)calloc(n,sizeof(double*));
  for (i = 0; i < n; i++)
    L[i] = (double*)calloc(n,sizeof(double));

  //now, init to zeros
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      L[i][j] = 0.0;         
  
  //now, fill in L with lower tri A, ones on diag
  for (i = 0; i < n; i++)
    {
    for (j = 0; j < n; j++)
      {
      if (i > j)
        L[i][j] = A[i][j];
      if (i == j)
        L[i][j] = 1.0;
      }
    }

  //now, do forward subs
  for (i = 0; i < n; i++)
    {
    x[i] /= L[i][i];
    for (j = i+1; j < n; j++)
      x[j] -= x[i]*L[j][i]; 
    }

  //clean up memory
  for (i = 0; i < n; i++)
    free(L[i]);
  freenull(L);

  return;
  }
