#include "Lu4.h"

void Lu4(double **A, const int n)
  {
  int i, j, k, m; //counters

  //first, take lower triangular elements and divide by diagonal elements in each column
  for (j = 0; j < n; j++)
    {
    for (i = j; i < n-1; i++)
      {
      A[i+1][j] /= A[j][j]; 
      }
    //now that we have done a column, we operate on columns and rows beyond that
    for (k = j+1; k < n; k++)
      {
      for (m = j+1; m < n; m++)
        {
        A[m][k] -= A[m][j]*A[j][k];
        }
      } 
    }
  return;
  }
