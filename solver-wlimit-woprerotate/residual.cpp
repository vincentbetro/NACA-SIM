#include "Mesh.h"
#include "analyticjacobian.h"
#include "splitflux.h"
//#include "time.h"

/*//for debugging
void wait ( int seconds )
  {
  clock_t endwait;
  endwait = clock () + seconds * CLOCKS_PER_SEC ;
  while (clock() < endwait) {}
  }*/

//we'll also do dR/dQ here and place in matrix to avoid redundant computations
void Mesh::residual()
  {
  int n0, n1, n2, nd0, nd1, nd2; //node placeholders to avoid extra memory access
  int i, j, k, l, m, flag; //counters
  Point corner1, corner2, corner3, corner4, p0, p1, p2, tricent, secmdpt; //point class placeholders, corners, centroid, second order cv side mdpt
  Vector u1, v1, norm, s0, s1, s2, norm0, norm1, norm2, norm0m, norm1m, norm2m, r1, r2; //vector class placeholders, side vectors, norms of side vectors, midpoint versions too, second order vectors
  Point mdpt0, mdpt1, mdpt2; //three edge midpoints for ease of use
  double l0, l1, l2, lm0, lm1, lm2, xnorm, ynorm, rlen, pl, pr, ufr, vfr, ufl, vfl, qtempL, qtempR, qtemp, phitemp; //edge lengths, normals for switch, left and right state dummy vars for bd
  
  //now, allocate for RHS matrix (actually a vector 4*nn, but for ease of access, we make 2-D)
  //based on previous q
  if (RHS == 0)
    {
    RHS = (double**)calloc(nn,sizeof(double*));
    for (i = 0; i < nn; i++)
      RHS[i] = (double*)calloc(4,sizeof(double));
    //printf("\nAllocating for RHS.");
    }
  else
    {
    for (i = 0; i < nn; i++)
      for (j = 0; j < 4; j++)
        RHS[i][j] = 0.0; //if already allocated, just reinit
    //printf("\nOnly redimensioning RHS.");
    }

  //allocate temp arrays for qr, ql, fp, fm, afp, afm to be passed in and out of split flux routine, analytic jacobian routine
  double *ql, *qr, *fps, *fms, **afp, **afm;
  ql = (double*)calloc(4,sizeof(double));
  qr = (double*)calloc(4,sizeof(double));
  fps = (double*)calloc(4,sizeof(double));
  fms = (double*)calloc(4,sizeof(double));
  if (impex == 1)
    {
    afp = (double**)calloc(4,sizeof(double*));
    afm = (double**)calloc(4,sizeof(double*));
    for (i = 0; i < 4; i++)
      {
      afp[i] = (double*)calloc(4,sizeof(double));
      afm[i] = (double*)calloc(4,sizeof(double));
      }
    }

  double *pl_ql, *pr_qr; //declare regardless for complier
  if (impex == 1)
    {
    //create temp arrays of p wrt ql and qr to facilitate flux jacobians for invscid bd
    pl_ql = (double*)calloc(4,sizeof(double));
    pr_qr = (double*)calloc(4,sizeof(double));
    } 

  //now, create temp arrays of gradient vectors
  Vector *grad1, *grad2;
  if (impex == 1 && order == 2)
    {
    grad1 = new Vector[4];
    grad2 = new Vector[4];
    } 

  flag = 0;
  //now, we loop thru the triangles and compute F+(QL) and F-(QR) and add them to get F on each face and add or subtract from correct node's RHS
  for (i = 0; i < nt && !flag; i++)
    {
    if (impex != 1 && order != 2)
      {
      flag = 1;
      continue;
      }
    
    //use node indices for cleanliness
    n0 = tri[i][0];
    n1 = tri[i][1];
    n2 = tri[i][2];
      
    //create point class placeholders for ease
    p0 = Point(nodep[n0][0],nodep[n0][1]);
    p1 = Point(nodep[n1][0],nodep[n1][1]);
    p2 = Point(nodep[n2][0],nodep[n2][1]);

    //create vectors from each side
    s0 = Vector(p0,p1);
    s1 = Vector(p1,p2);
    s2 = Vector(p2,p0);
    //find the sides' magnitudes
    l0 = s0.magnitude();
    l1 = s1.magnitude();
    l2 = s2.magnitude();
    //find normals
    norm0 = Vector(s0[1],-s0[0]);
    norm1 = Vector(s1[1],-s1[0]);
    norm2 = Vector(s2[1],-s2[0]);
    //normalize them
    norm0 /= l0;
    norm1 /= l1;
    norm2 /= l2;

    //find triangle center by averaging the nodes
    tricent = (p0 + p1 + p2)/3.0;
    //printf("\n(%lf, %lf)\n",tricent[0],tricent[1]);

    //find midpoints
    mdpt0 = (p0 + p1)/2.0;
    mdpt1 = (p1 + p2)/2.0;
    mdpt2 = (p2 + p0)/2.0;

    //create vectors from each cv piece
    s0 = Vector(mdpt0,tricent);
    s1 = Vector(mdpt1,tricent);
    s2 = Vector(mdpt2,tricent);
    //find the cv pieces' magnitudes
    lm0 = s0.magnitude();
    lm1 = s1.magnitude();
    lm2 = s2.magnitude();
    //find normals
    norm0m = Vector(s0[1],-s0[0]);
    norm1m = Vector(s1[1],-s1[0]);
    norm2m = Vector(s2[1],-s2[0]);
    //normalize them
    norm0m /= lm0;
    norm1m /= lm1;
    norm2m /= lm2;

    //now, find the split fluxes for each node using ql and qr, sum them into flux, and begin to sum/subtract them as appropriate to RHS
    //we use a switch
    for (j = 0; j < 3; j++)
      {
      switch(j)
        {
        case 0:
            secmdpt = (mdpt0 + tricent)/2.0;
            r1 = Vector(p0,secmdpt);
            r2 = Vector(p1,secmdpt);       
            //first, determine ql, qr
            for (k = 0; k < 4; k++)
              {
              grad1[k] = Vector(dqdx[n0][k],dqdy[n0][k]);
              grad2[k] = Vector(dqdx[n1][k],dqdy[n1][k]);
              
              //temp vars
              qtempL = q[n0][k] + grad1[k]*r1;
              qtempR = q[n1][k] + grad2[k]*r2;
              
              //check for limiting..if == , phi init to 1.0
              if (qtempL > q[n0][k])
                {
                qtemp = (qmax[n0][k] - q[n0][k])/(qtempL - q[n0][k]);
                phitemp = MIN(1.0,qtemp);
                phi[n0][k] = MIN(phi[n0][k],phitemp);
                }
              
              if (qtempL < q[n0][k])
                {
                qtemp = (qmin[n0][k] - q[n0][k])/(qtempL - q[n0][k]);
                phitemp = MIN(1.0,qtemp);
                phi[n0][k] = MIN(phi[n0][k],phitemp);
                }
                
              if (qtempR > q[n1][k])
                {
                qtemp = (qmax[n1][k] - q[n1][k])/(qtempR - q[n1][k]);
                phitemp = MIN(1.0,qtemp);
                phi[n1][k] = MIN(phi[n1][k],phitemp);
                }
              
              if (qtempR < q[n1][k])
                {
                qtemp = (qmin[n1][k] - q[n1][k])/(qtempR - q[n1][k]);
                phitemp = MIN(1.0,qtemp);
                phi[n1][k] = MIN(phi[n1][k],phitemp);
                }
              }     
        break;
        case 1:
            secmdpt = (mdpt1 + tricent)/2.0;
            r1 = Vector(p1,secmdpt);
            r2 = Vector(p2,secmdpt);
            //first, determine ql, qr
            for (k = 0; k < 4; k++)
              {
              grad1[k] = Vector(dqdx[n1][k],dqdy[n1][k]);
              grad2[k] = Vector(dqdx[n2][k],dqdy[n2][k]);
              
              //temp vars
              qtempL = q[n1][k] + grad1[k]*r1;
              qtempR = q[n2][k] + grad2[k]*r2;
              
              //check for limiting..if == , phi init to 1.0
              if (qtempL > q[n1][k])
                {
                qtemp = (qmax[n1][k] - q[n1][k])/(qtempL - q[n1][k]);
                phitemp = MIN(1.0,qtemp);
                phi[n1][k] = MIN(phi[n1][k],phitemp);
                }
              
              if (qtempL < q[n1][k])
                {
                qtemp = (qmin[n1][k] - q[n1][k])/(qtempL - q[n1][k]);
                phitemp = MIN(1.0,qtemp);
                phi[n1][k] = MIN(phi[n1][k],phitemp);
                }
                
              if (qtempR > q[n2][k])
                {
                qtemp = (qmax[n2][k] - q[n2][k])/(qtempR - q[n2][k]);
                phitemp = MIN(1.0,qtemp);
                phi[n2][k] = MIN(phi[n2][k],phitemp);
                }
              
              if (qtempR < q[n2][k])
                {
                qtemp = (qmin[n2][k] - q[n2][k])/(qtempR - q[n2][k]);
                phitemp = MIN(1.0,qtemp);
                phi[n2][k] = MIN(phi[n2][k],phitemp);
                }
              }     
        break;
        case 2:
            secmdpt = (mdpt2 + tricent)/2.0;
            r1 = Vector(p2,secmdpt);
            r2 = Vector(p0,secmdpt);
            //first, determine ql, qr
            for (k = 0; k < 4; k++)
              {
              grad1[k] = Vector(dqdx[n2][k],dqdy[n2][k]);
              grad2[k] = Vector(dqdx[n0][k],dqdy[n0][k]);
              
              //temp vars
              qtempL = q[n2][k] + grad1[k]*r1;
              qtempR = q[n0][k] + grad2[k]*r2;
              
              //check for limiting..if == , phi init to 1.0
              if (qtempL > q[n2][k])
                {
                qtemp = (qmax[n2][k] - q[n2][k])/(qtempL - q[n2][k]);
                phitemp = MIN(1.0,qtemp);
                phi[n2][k] = MIN(phi[n2][k],phitemp);
                }
              
              if (qtempL < q[n2][k])
                {
                qtemp = (qmin[n2][k] - q[n2][k])/(qtempL - q[n2][k]);
                phitemp = MIN(1.0,qtemp);
                phi[n2][k] = MIN(phi[n2][k],phitemp);
                }
                
              if (qtempR > q[n0][k])
                {
                qtemp = (qmax[n0][k] - q[n0][k])/(qtempR - q[n0][k]);
                phitemp = MIN(1.0,qtemp);
                phi[n0][k] = MIN(phi[n0][k],phitemp);
                }
              
              if (qtempR < q[n0][k])
                {
                qtemp = (qmin[n0][k] - q[n0][k])/(qtempR - q[n0][k]);
                phitemp = MIN(1.0,qtemp);
                phi[n0][k] = MIN(phi[n0][k],phitemp);
                }
              }
        break;
        default:
          printf("\nYou have a triangle with more than three sides according to the loop counter in residual.cpp!\n");
          fflush(stdout);
          exit(0);    
        break;
        }
      }
    }
    
    
    for (i = 0; i < nt; i++)
    {
    //use node indices for cleanliness
    n0 = tri[i][0];
    n1 = tri[i][1];
    n2 = tri[i][2];
      
    //create point class placeholders for ease
    p0 = Point(nodep[n0][0],nodep[n0][1]);
    p1 = Point(nodep[n1][0],nodep[n1][1]);
    p2 = Point(nodep[n2][0],nodep[n2][1]);

    //create vectors from each side
    s0 = Vector(p0,p1);
    s1 = Vector(p1,p2);
    s2 = Vector(p2,p0);
    //find the sides' magnitudes
    l0 = s0.magnitude();
    l1 = s1.magnitude();
    l2 = s2.magnitude();
    //find normals
    norm0 = Vector(s0[1],-s0[0]);
    norm1 = Vector(s1[1],-s1[0]);
    norm2 = Vector(s2[1],-s2[0]);
    //normalize them
    norm0 /= l0;
    norm1 /= l1;
    norm2 /= l2;

    //find triangle center by averaging the nodes
    tricent = (p0 + p1 + p2)/3.0;
    //printf("\n(%lf, %lf)\n",tricent[0],tricent[1]);

    //find midpoints
    mdpt0 = (p0 + p1)/2.0;
    mdpt1 = (p1 + p2)/2.0;
    mdpt2 = (p2 + p0)/2.0;

    //create vectors from each cv piece
    s0 = Vector(mdpt0,tricent);
    s1 = Vector(mdpt1,tricent);
    s2 = Vector(mdpt2,tricent);
    //find the cv pieces' magnitudes
    lm0 = s0.magnitude();
    lm1 = s1.magnitude();
    lm2 = s2.magnitude();
    //find normals
    norm0m = Vector(s0[1],-s0[0]);
    norm1m = Vector(s1[1],-s1[0]);
    norm2m = Vector(s2[1],-s2[0]);
    //normalize them
    norm0m /= lm0;
    norm1m /= lm1;
    norm2m /= lm2;

    //now, find the split fluxes for each node using ql and qr, sum them into flux, and begin to sum/subtract them as appropriate to RHS
    //we use a switch
    for (j = 0; j < 3; j++)
      {
      switch(j)
        {
        case 0:
          if (impex == 1 && order == 2)
            {
            secmdpt = (mdpt0 + tricent)/2.0;
            r1 = Vector(p0,secmdpt);
            r2 = Vector(p1,secmdpt);       
            //first, determine ql, qr
            for (k = 0; k < 4; k++)
              {
              grad1[k] = Vector(dqdx[n0][k],dqdy[n0][k]);
              grad2[k] = Vector(dqdx[n1][k],dqdy[n1][k]);
              
              ql[k] = q[n0][k] + phi[n0][k]*(grad1[k]*r1);
              qr[k] = q[n1][k] + phi[n1][k]*(grad2[k]*r2);
              }
            //determine xnorm and ynorm, rlen
            xnorm = norm0m[0];
            ynorm = norm0m[1];
            rlen = lm0;
            //set node placeholders
            nd0 = n0;
            nd1 = n1;
            }
          else
            {     
            //first, determine ql, qr
            for (k = 0; k < 4; k++)
              {
              ql[k] = q[n0][k];
              qr[k] = q[n1][k];
              }
            //determine xnorm and ynorm, rlen
            xnorm = norm0m[0];
            ynorm = norm0m[1];
            rlen = lm0;
            //set node placeholders
            nd0 = n0;
            nd1 = n1;
            }     
        break;
        case 1:
          if (impex == 1 && order == 2)
            {
            secmdpt = (mdpt1 + tricent)/2.0;
            r1 = Vector(p1,secmdpt);
            r2 = Vector(p2,secmdpt);
            //first, determine ql, qr
            for (k = 0; k < 4; k++)
              {
              grad1[k] = Vector(dqdx[n1][k],dqdy[n1][k]);
              grad2[k] = Vector(dqdx[n2][k],dqdy[n2][k]);
              
              ql[k] = q[n1][k] + phi[n1][k]*(grad1[k]*r1);
              qr[k] = q[n2][k] + phi[n2][k]*(grad2[k]*r2);
              }
            //determine xnorm and ynorm, rlen
            xnorm = norm1m[0];
            ynorm = norm1m[1];
            rlen = lm1;
            //set node placeholders
            nd0 = n1;
            nd1 = n2;
            }
          else
            {     
            //first, determine ql, qr
            for (k = 0; k < 4; k++)
              {
              ql[k] = q[n1][k];
              qr[k] = q[n2][k];
              }
            //determine xnorm and ynorm, rlen
            xnorm = norm1m[0];
            ynorm = norm1m[1];
            rlen = lm1;
            //set node placeholders
            nd0 = n1;
            nd1 = n2;
            }     
        break;
        case 2:
          if (impex == 1 && order == 2)
            {
            secmdpt = (mdpt2 + tricent)/2.0;
            r1 = Vector(p2,secmdpt);
            r2 = Vector(p0,secmdpt);
            //first, determine ql, qr
            for (k = 0; k < 4; k++)
              {
              grad1[k] = Vector(dqdx[n2][k],dqdy[n2][k]);
              grad2[k] = Vector(dqdx[n0][k],dqdy[n0][k]);
              
              ql[k] = q[n2][k] + phi[n2][k]*(grad1[k]*r1);
              qr[k] = q[n0][k] + phi[n0][k]*(grad2[k]*r2);
              }
            //determine xnorm and ynorm, rlen
            xnorm = norm2m[0];
            ynorm = norm2m[1];
            rlen = lm2;
            //set node placeholders
            nd0 = n2;
            nd1 = n0;
            }
          else
            {     
            //first, determine ql, qr
            for (k = 0; k < 4; k++)
              {
              ql[k] = q[n2][k];
              qr[k] = q[n0][k];
              }
            //determine xnorm and ynorm, rlen
            xnorm = norm2m[0];
            ynorm = norm2m[1];
            rlen = lm2;
            //set node placeholders
            nd0 = n2;
            nd1 = n0;
            }     
        break;
        default:
          printf("\nYou have a triangle with more than three sides according to the loop counter in residual.cpp!\n");
          fflush(stdout);
          exit(0);    
        break;
        }
      //now, determine the split fluxes
      splitflux(ql, qr, gamma, xnorm, ynorm, rlen, fps, fms);
      if (impex == 1)
        {
        //now, determine analytic jacobians
        analyticjacobian(ql, qr, gamma, xnorm, ynorm, rlen, fps, fms, afp, afm);
        }
      //now, sum them to get F (hold result in fp to save space)
      for (k = 0; k < 4; k++)
        {
        fps[k] += fms[k];
        }
      //now, place in RHS approproiately (remember, we are subtracting over, so signs are opposite)
      for (k = 0; k < 4; k++)
        {
        RHS[nd0][k] -= fps[k];
        RHS[nd1][k] += fps[k];
        }
      if (impex == 1)
        {
        //now, place afp and afm in proper diagonal entries in mat
        for (k = 0; k < 4; k++)
          {
          for (l = 0; l < 4; l++)
            {
            //if free stream, do not add to diagonal
            if ((tag[nd0] == 0 && tagbd[nd0] == 1) || (tag[nd0] == 1 && tagbd[nd0] == 0))
              mat[iau[nd0]][k][l] += afp[k][l];
            //if free stream, do not add to diagonal
            if ((tag[nd1] == 0 && tagbd[nd1] == 1) || (tag[nd1] == 1 && tagbd[nd1] == 0))
              mat[iau[nd1]][k][l] -= afm[k][l];         
            }
          }
        //now, place afp and afm in proper off diagonal entries in mat
        //if free stream, do not add to off diagonal
        if ((tag[nd0] == 0 && tagbd[nd0] == 1) || (tag[nd0] == 1 && tagbd[nd0] == 0))
          {
          for (m = ia[nd0]; m < ia[nd0+1]; m++)
            {
            if (ja[m] != nd1)
              continue; //this is the diagonal element
            for (k = 0; k < 4; k++)
              {
              for (l = 0; l < 4; l++)
                {
                mat[m][k][l] += afm[k][l];
                }
              }
            }
          }
        //if free stream, do not add to off diagonal
        if ((tag[nd1] == 0 && tagbd[nd1] == 1) || (tag[nd1] == 1 && tagbd[nd1] == 0))
          {
          for (m = ia[nd1]; m < ia[nd1+1]; m++)
            {
            if (ja[m] != nd0)
              continue; //this is the diagonal element
            for (k = 0; k < 4; k++)
              {
              for (l = 0; l < 4; l++)
                {
                mat[m][k][l] -= afp[k][l];
                }
              }
            }
          }
        } //end implicit info
      }

    //now, we have to add in boundary contributions (inner, inviscid bd only, since outer bd we have flux (and RHS) = 0)
    for (j = 0; j < 3; j++)
      {
      switch(j)
        {
        case 0:
          s0 = Vector(p0,mdpt0);
          s1 = Vector(mdpt0,p1);
          l0 = s0.magnitude();
          l1 = s1.magnitude();
          norm0 = Vector(s0[1],-s0[0]);
          norm1 = Vector(s1[1],-s1[0]);
          norm0 /= l0;
          norm1 /= l1;
          nd0 = n0;
          nd1 = n1;  
        break;
        case 1:
          s0 = Vector(p1,mdpt1);
          s1 = Vector(mdpt1,p2);
          l0 = s0.magnitude();
          l1 = s1.magnitude();
          norm0 = Vector(s0[1],-s0[0]);
          norm1 = Vector(s1[1],-s1[0]);
          norm0 /= l0;
          norm1 /= l1;
          nd0 = n1;
          nd1 = n2; 
        break;
        case 2:
          s0 = Vector(p2,mdpt2);
          s1 = Vector(mdpt2,p0);
          l0 = s0.magnitude();
          l1 = s1.magnitude();
          norm0 = Vector(s0[1],-s0[0]);
          norm1 = Vector(s1[1],-s1[0]);
          norm0 /= l0;
          norm1 /= l1;
          nd0 = n2;
          nd1 = n0; 
        break;
        default:
          printf("\nYou have a triangle with more than three sides according to the loop counter in residual.cpp!\n");
          fflush(stdout);
          exit(0);
        break;
        }
       
     //check if we have a boundary node, do calculations on given side
     //note that first and last entries of flux are 0, so we only manipulate RHS [1], [2]
     if (tagbd[nd0] == 1 && tagbd[nd1] == 1)
       {
       //to make things easier to read
       pl = (gamma - 1.0)*(q[nd0][3] - 0.5*((q[nd0][1]*q[nd0][1] + q[nd0][2]*q[nd0][2])/q[nd0][0]));
       pr = (gamma - 1.0)*(q[nd1][3] - 0.5*((q[nd1][1]*q[nd1][1] + q[nd1][2]*q[nd1][2])/q[nd1][0]));

       RHS[nd0][1] -= 0.75*(norm0[0]*l0*pl) + 0.25*(norm1[0]*l1*pr);
       RHS[nd0][2] -= 0.75*(norm0[1]*l0*pl) + 0.25*(norm1[1]*l1*pr);
       RHS[nd1][1] -= 0.75*(norm1[0]*l1*pr) + 0.25*(norm0[0]*l0*pl);
       RHS[nd1][2] -= 0.75*(norm1[1]*l1*pr) + 0.25*(norm0[1]*l0*pl);

       if (impex == 1)
         {
         //now, we compute afp and afm for special flux
         //1st, 4th rows are all zeros
         for (k = 0; k < 4; k++)
           afp[0][k] = afp[3][k] = afm[0][k] = afm[3][k] = 0.0;

         //for creating diff p, use weighted values from both nodes
         ufl = 0.75*(q[nd0][1]/q[nd0][0]) + 0.25*(q[nd1][1]/q[nd1][0]);
         vfl = 0.75*(q[nd0][2]/q[nd0][0]) + 0.25*(q[nd1][2]/q[nd1][0]);

         ufr = 0.75*(q[nd1][1]/q[nd1][0]) + 0.25*(q[nd0][1]/q[nd0][0]);
         vfr = 0.75*(q[nd1][2]/q[nd1][0]) + 0.25*(q[nd0][2]/q[nd0][0]); 
       
         //set up differentiated p wrt q0..3
         pl_ql[0] = 0.5*(gamma - 1.0)*(ufl*ufl + vfl*vfl);
         pl_ql[1] = (1.0 - gamma)*ufl;
         pl_ql[2] = (1.0 - gamma)*vfl;
         pl_ql[3] = (gamma - 1.0);

         pr_qr[0] = 0.5*(gamma - 1.0)*(ufr*ufr + vfr*vfr);
         pr_qr[1] = (1.0 - gamma)*ufr;
         pr_qr[2] = (1.0 - gamma)*vfr;
         pr_qr[3] = (gamma - 1.0);
    
         //compute flux jacobians
         for (k = 0; k < 4; k++)
           {
           afp[1][k] = norm0[0]*l0*pl_ql[k];
           afp[2][k] = norm0[1]*l0*pl_ql[k];
           afm[1][k] = norm1[0]*l1*pr_qr[k];
           afm[2][k] = norm1[1]*l1*pr_qr[k];
           }
      
         //now, place afp and afm in proper diagonal entries in mat
         for (k = 0; k < 4; k++)
           {
           for (l = 0; l < 4; l++)
             {
             mat[iau[nd0]][k][l] += afp[k][l];
             mat[iau[nd1]][k][l] += afm[k][l];         
             }
           }
         } //end implicit info
       }
     }
   }
  //now, after all is said and done, reset RHS for far field to zero    
  for (i = 0; i < nt; i++)
    {
    //use node indices for cleanliness
    n0 = tri[i][0];
    n1 = tri[i][1];
    n2 = tri[i][2];
    
    //switch thru sides
    for (j = 0; j < 3; j++)
      {
      switch(j)
        {
        case 0:
          //set node placeholders
          nd0 = n0;
          nd1 = n1;     
        break;
        case 1:
          //set node placeholders
          nd0 = n1;
          nd1 = n2;     
        break;
        case 2:
          //set node placeholders
          nd0 = n2;
          nd1 = n0;     
        break;
        default:
          printf("\nYou have a triangle with more than three sides according to the loop counter in residual.cpp!\n");
          fflush(stdout);
          exit(0);    
        break;
        }
      //actually reset RHS
      if (tag[nd0] == 0 && tag[nd1] == 0 && tagbd[nd0] == 0 && tagbd[nd1] == 0)
        {
        for (k = 0; k < 4; k++)
          {
          RHS[nd0][k] = 0.0;
          RHS[nd1][k] = 0.0;
          }
        }
      }
    }


  /*//debug, create list
  List stops;
  int endlist = nn/10; //integer division
  for (i = 0; i < endlist; i++)
    stops.Check_List(i*10);
  
  for (i = 0; i < nn; i++)
    {
    if (stops.Is_In_List(i))
         wait(10);
    for (j = 0; j < 4; j++)
      {
      printf("RHS[%d][%d] = %lf\n",i,j,RHS[i][j]);
      }
    }*/
       
  //clean up memory
  freenull(ql);
  freenull(qr);
  freenull(fps);
  freenull(fms);
  if (impex == 1)
    {
    for (i = 0; i < 4; i++)
      {
      free(afp[i]);
      free(afm[i]);
      } 
    freenull(afp);
    freenull(afm);
    freenull(pl_ql);
    freenull(pr_qr);
    }
  if (impex == 1 && order == 2)
    {
    delete [] grad1;
    delete [] grad2;
    }

  return;
  }
