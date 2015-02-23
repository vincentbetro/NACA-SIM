#include "Mesh.h"

void Mesh::mat_build()
  {
  int i, j, k; //counters

  double deltaT; //time step for each node will vary

  double p, magubar, rho; //broken out derived vars for ease

  //now, allocate for mat (holds 4x4 matrices for each node, reinit to zero each time)
  if (mat == 0)
    {
    mat = (double***)calloc(mdim,sizeof(double**));
    for (i = 0; i < mdim; i++)
      {
      mat[i] = (double**)calloc(4,sizeof(double*));
      for (j = 0; j < 4; j++)
        {
        mat[i][j] = (double*)calloc(4,sizeof(double));
        }
      }
    //printf("\nAllocating for mat.");
    }
  else
    {
    for (i = 0; i < mdim; i++)
      for (j = 0; j < 4; j++)
        for (k = 0; k < 4; k++)
          mat[i][j][k] = 0.0; //just reinit, but do not realloc
    //printf("\nOnly redimensioning mat.");
    }
   
  //printf("\nCFL this run (mesh_obj) = %lf\n",CFL);
      
  //now, we build the global matrix..for now, only diagonal entries of diagonal matrices are filled in, the dR/dQ filled in in residual.cpp to avoid excess comp
  for (i = 0; i < nn; i++)
    {
    //for far field
    if (tag[i] == 0 && tagbd[i] == 0)
      {
      mat[iau[i]][0][0] = 1.0;
      mat[iau[i]][0][1] = 0.0;
      mat[iau[i]][0][2] = 0.0;
      mat[iau[i]][0][3] = 0.0;
      mat[iau[i]][1][0] = 0.0;
      mat[iau[i]][1][1] = 1.0;
      mat[iau[i]][1][2] = 0.0;
      mat[iau[i]][1][3] = 0.0;
      mat[iau[i]][2][0] = 0.0;
      mat[iau[i]][2][1] = 0.0;
      mat[iau[i]][2][2] = 1.0;
      mat[iau[i]][2][3] = 0.0;
      mat[iau[i]][3][0] = 0.0;
      mat[iau[i]][3][1] = 0.0;
      mat[iau[i]][3][2] = 0.0;
      mat[iau[i]][3][3] = 1.0;
      }
    else
      {
      //for inviscid and internal
      //first, compute deltaT
      //break out p, mag(ubar), rho for ease
      p = (gamma - 1.0)*(q[i][3] - 0.5*((q[i][1]*q[i][1] + q[i][2]*q[i][2])/q[i][0]));
      magubar = sqrt((q[i][1]/q[i][0])*(q[i][1]/q[i][0]) + (q[i][2]/q[i][0])*(q[i][2]/q[i][0]));
      rho = q[i][0];

      deltaT = (CFL*len_scale[i])/(magubar + sqrt((gamma*p)/rho));

      if (fabs(deltaT) > 1.0e5 || fabs(deltaT) < 1.0e-15)
        printf("\ndeltaT[%d] with tag = %d is bad.\n",i,tag[i]);

      if (fabs(cvareaG[i]) > 1.0e5 || fabs(cvareaG[i]) < 1.0e-15)
        printf("\ncvareaG[%d] with tag = %d is bad.\n",i,tag[i]);

      mat[iau[i]][0][0] = cvareaG[i]/deltaT;
      mat[iau[i]][0][1] = 0.0;
      mat[iau[i]][0][2] = 0.0;
      mat[iau[i]][0][3] = 0.0;
      mat[iau[i]][1][0] = 0.0;
      mat[iau[i]][1][1] = cvareaG[i]/deltaT;
      mat[iau[i]][1][2] = 0.0;
      mat[iau[i]][1][3] = 0.0;
      mat[iau[i]][2][0] = 0.0;
      mat[iau[i]][2][1] = 0.0;
      mat[iau[i]][2][2] = cvareaG[i]/deltaT;
      mat[iau[i]][2][3] = 0.0;
      mat[iau[i]][3][0] = 0.0;
      mat[iau[i]][3][1] = 0.0;
      mat[iau[i]][3][2] = 0.0;
      mat[iau[i]][3][3] = cvareaG[i]/deltaT;
      }
    }

   
  return;
  } 
 
