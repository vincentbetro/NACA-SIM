#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Point.h"
#include "Vector.h"
#include "List.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
#define freenull(a) {if(a) free(a); a=NULL;}

#ifndef Mesh_h
#define Mesh_h
class Mesh
  {
  public:
    Mesh()
      { 
      //init vars in constructor
      mdim = nn = nt = nq = nb = nc = nv = 0; //dim of mat, num nodes, tris, quads, bd, num constants, num vars
      nblk = 1; //num blocks
      //limiter
      phi = 0;
      qmax = 0;
      qmin = 0;
      //arrays of values read in
      nodep = 0;
      tri = 0;
      quad = 0;
      constants = 0;
      //boundary arrays
      nbs = 0;
      bs = 0;
      //hash table
      hash = 0;
      //compressed row storage arrays
      ia = 0;
      ja = 0;
      iau = 0;
      //A matrix
      mat = 0;
      //RHS (residual)
      RHS = 0;
      //control volume area using green's theorem, length scale
      cvareaG = 0;
      len_scale = 0;
      //q at nodes
      q = 0;
      //grad at nodes
      dqdx = 0;
      dqdy = 0;
      //bd nodes tag array
      tag = 0;
      //inner or outer bd tag array
      tagbd = 0;
      //relaxation factor init
      omega = 1.0;
      //cvg value init
      cvg = 1.0e-6;
      //gamma init
      gamma = 1.4;
      //CFL init
      CFL = 1.0;
      //inner boundary init
      ib1 = ib2 = 0;
      //init implicit/explicit toggle, order toggle
      impex = 0;
      order = 0;
      }
    ~Mesh() 
      {
      int i, j;
      //relaxation factor re-init
      omega = 0.0;
      //cvg value re-init
      cvg = 0.0;
      //gamma re-init
      gamma = 0.0;
      //CFL re-init
      CFL = 0.0;
      //inner boundary re-init
      ib1 = ib2 = 0;
      //limiter
      if (phi != 0)
        {
        for (i = 0; i < nn; i++)
          free(phi[i]);
        freenull(phi);
        }
      if (qmax != 0)
        {
        for (i = 0; i < nn; i++)
          free(qmax[i]);
        freenull(qmax);
        }
      if (qmin != 0)
        {
        for (i = 0; i < nn; i++)
          free(qmin[i]);
        freenull(qmin);
        }
      //arrays of values read in
      for (i = 0; i < nn; i++)
        free(nodep[i]);
      freenull(nodep);
      for (i = 0; i < nt; i++)
        free(tri[i]);
      freenull(tri);
      if (nq > 0)
        {
        for (i = 0; i < nq; i++)
          free(quad[i]);
        freenull(quad);
        }
      if (nc > 0)
        freenull(constants);
      //boundary arrays
      for (i = 0; i < nb; i++)
        for (j = 0; j < nbs[i]; j++)
          free(bs[i][j]);
      for (i = 0; i < nbs[i]; i++)
        free(bs[i]);
      freenull(bs);
      freenull(nbs);
      //hash table
      if (hash != 0)
        {
        for (i = 0; i < nn; i++)
          delete hash[i];
        delete [] hash;
        }
      //compressed row storage arrays
      freenull(ia);
      freenull(ja);
      freenull(iau);
      //A matrix
      for (i = 0; i < mdim; i++)
        for (j = 0; j < 4; j++)
          free(mat[i][j]);
      for (i = 0; i < 4; i++)
        free(mat[i]);
      freenull(mat);
      //RHS (residual)
      for (i = 0; i < nn; i++)
        free(RHS[i]);
      freenull(RHS);
      //control volume area using green's theorem, length scale
      freenull(cvareaG);
      freenull(len_scale);
      //q at nodes
      for (i = 0; i < nn; i++)
        free(q[i]);
      freenull(q);
      //grad at nodes
      if (impex == 1 && order == 2)
        { 
        for (i = 0; i < nn; i++)
          {
          free(dqdx[i]);
          free(dqdy[i]);
          }
        freenull(dqdx);
        freenull(dqdy);
        }
      //bd nodes tag array
      freenull(tag);
      //inner or outer bd tag array
      freenull(tagbd);
      //finally, after use in destructor, re-init vars in destructor (to make obvious that it destructed)
      mdim = nblk = nn = nt = nq = nb = nc = nv = 0; //dim of mat, num blocks, num nodes, tris, quads, bd
      //reset implicit/explicit toggle, order toggle
      impex = 0;
      order = 0;
      }

  //declare functions
  void mesh_read(char filename[]); //reads in mesh
  void mesh_out(char filename[]); //writes out mesh
  void mesh_init(double rho0, double u0, double v0, double et0); //computes areas with green's theorem, inits q
  void Lhash(); //creates hash table, creates mdim
  void cuthill_mckee(); //does cuthill-mckee reordering to make implicit methods faster by bunching points along diagonal
  void crs(); //creates compressed row storage
  void residual(); //creates RHS vector
  void gradients(); //computes gradients with green's theorem
  void mat_build(); //builds A matrix (mat)
  double lin_solve(int &tempiter); //takes max internal iterations, passes out acutal, returns RMS error

  //delcare vars
  int mdim, nn, nblk, nt, nq, nb, ib1, ib2, nc, nv; //dim of mat, number of nodes, number of blocks, number of triangles, number of quads, number of boundaries, inner boundaries, num of const, num vars
  int impex, order; //toggles for types of solve
  double omega, cvg, gamma, CFL;  //read in values 
  //declare matrices for x,y coords of physical nodes (indexed), triangle conn, quad conn
  double **nodep;
  int **tri, **quad;
  double *constants; //constants read in from mesh file 
  //boundary segments, indexed as boundary number, edge number, node number on edge; number boundary segments per edge 
  int ***bs, *nbs;
  //areas of control volumes as well as q and gradients of q at each node, length scale
  double *cvareaG, **q, **dqdx, **dqdy, *len_scale;
  //tag array of boundary nodes, inner or outer is tagbd
  int *tag, *tagbd;
  //hash table
  List **hash;
  //compressed row storage arrays
  int *ia, *ja, *iau;
  //A matrix
  double ***mat;
  //RHS (residual)
  double **RHS;
  //limiter 
  double **phi, **qmax, **qmin;

  };

#endif
