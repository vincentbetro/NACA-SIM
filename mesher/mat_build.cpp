#include "mat_build.h"

//we pass in ux, uy, vx, vy to compute winslow coeff (these are sent in as 1 element arrays from L-E and ignored), poisson for L-E (set to zero for winslow), phi and psi (forcing functions), number of nodes, reinint to zeros version of mat to come back full, number of triangles, tri, from winslow, pass in nodec for node, from linear_elastic, pass in nodep for node, and compressed row storage indexing arrays 
//the rub that won't let us compute these externally is that the alphas and betas are node based in winslow and the alphas and thetas are element based in L-E...this is why we must pass in smooth_type 
//pass in weights for forcing functions and overall forcing function weight (C)
void mat_build(int smooth_type, double *ux, double *uy, double *vx, double *vy, double poisson, double *phi, double *psi, int nn, double ***mat, int nt, int **tri, double **node, int *ia, int *ja, int *iau, double *w1, double *w2, double *w1x, double *w1y, double *w2x, double *w2y, double C)
  {

  int i, j, k; //counters
  int n0, n1, n2, nd0, nd1, nd2; //node variables to keep things clean and save on memory calls
  Vector u1, v1, norm0, norm1, norm2, s0, s1, s2, norm0a, norm1a, norm2a; //vectors to compute triangle area, normals, sides (coresponding to opposite node), normals for switch
  Point p0, p1, p2; //points to compose vectors
  double triarea, tx, ty; //area of given triangle when building weights, appropriate normals for coupling
  double w101, w102, w111, w112, w121, w122, w201, w202, w211, w212, w221, w222; //placeholders for weights..the first 1 or 2 stands for u or v, the second for the node index in the switch, and the third for the diag (1) or off-diag (2) in the submatrix
  double a, b, c, d, det; //placeholder vars for inverting diagonal matrices, determinant
  double alpha11, alpha12, alpha21, alpha22, theta11, theta12, theta21, theta22, beta; //winslow and L-E coefficients
  double E; //young's modulus
         
  //now, we build the global matrix
  for (i = 0; i < nt; i++)
    {
    //compute area first, since orientation doesn't change it (this is only duped since previous area calc was summative to nodes)
    //use initial n0, n1, n2 orientation for area calc, although it will be changed in switch
    n0 = tri[i][0];
    n1 = tri[i][1];
    n2 = tri[i][2];
    //create points, vectors to do area calulation
    p0 = Point(node[n0][0],node[n0][1]);
    p1 = Point(node[n1][0],node[n1][1]);
    p2 = Point(node[n2][0],node[n2][1]);
    u1 = Vector(p1,p2);
    v1 = Vector(p1,p0);
    //compute triarea
    triarea = .5*(u1 % v1);

    //now, re-compute norms 
    //create vectors of triangle sides (opposite node that is norm#)
    s0 = Vector(p1,p2);
    s1 = Vector(p2,p0);
    s2 = Vector(p0,p1);
    //create normal vectors
    norm0 = Vector(s0[1],-s0[0]);
    norm1 = Vector(s1[1],-s1[0]);
    norm2 = Vector(s2[1],-s2[0]);

    //if we are doing L-E, our coeff are element based, so we can compute outside j loop
    if (smooth_type == 1)
      {
      //here we must compute young_mod over the element in question
      E = young_mod(p0, p1, p2);

      //now, use this to compute weights
      alpha11 = alpha22 = (E*(1.0 - poisson))/((1.0 + poisson)*(1.0 - 2.0*poisson));
      alpha12 = alpha21 = E/(2.0*(1.0 + poisson));
      theta11 = theta22 = (E*poisson)/((1.0 + poisson)*(1.0 - 2.0*poisson));
      theta12 = theta21 = E/(2.0*(1.0 + poisson));
      beta = 0.0;  
      }

    // use a switch to work through nodes, this way, we always look from n0, but the definition changes to each node in triangle
    for (j = 0; j < 3; j++)
      {
      switch (j)
        {
        case 0: 
          //set the nodes on the triangle we are looking at from the perspective of each node
          nd0 = n0;
          nd1 = n1;
          nd2 = n2;
          tx = norm0[0]; //appropriate x-norm
          ty = norm0[1]; //appropriate y-norm
          //set normals based on which node we are looking at
          norm0a = norm0;
          norm1a = norm1;
          norm2a = norm2;
        break;
          
        case 1: 
          nd0 = n1;
          nd1 = n2;
          nd2 = n0;
          tx = norm1[0]; //appropriate x-norm
          ty = norm1[1]; //appropriate y-norm
          //set normals based on which node we are looking at
          norm0a = norm1;
          norm1a = norm2;
          norm2a = norm0; 
        break;
             
        case 2: 
          nd0 = n2;
          nd1 = n0;
          nd2 = n1;
          tx = norm2[0]; //appropriate x-norm
          ty = norm2[1]; //appropriate y-norm
          //set normals based on which node we are looking at
          norm0a = norm2;
          norm1a = norm0;
          norm2a = norm1;
        break;

        default: 
          printf("\nYour triangle is attempting to have more than 3 nodes in the switch in mat_build.cpp\n");
          fflush(stdout);
          exit(0); 
        break;
        }

      //we need to compute our coefficients for winslow, since they are node based, within j loop
      if (smooth_type == 0)
        {
        alpha11 = alpha21 = SIGN(MAX(fabs(uy[nd0]*uy[nd0] + vy[nd0]*vy[nd0]),1.0e-15),uy[nd0]*uy[nd0] + vy[nd0]*vy[nd0]);
        alpha12 = alpha22 = SIGN(MAX(fabs(ux[nd0]*ux[nd0] + vx[nd0]*vx[nd0]),1.0e-15),ux[nd0]*ux[nd0] + vx[nd0]*vx[nd0]);
        beta = SIGN(MAX(fabs(ux[nd0]*uy[nd0] + vx[nd0]*vy[nd0]),1.0e-15),ux[nd0]*uy[nd0] + vx[nd0]*vy[nd0]);
        theta11 = theta12 = theta21 = theta22 = 0.0;
        
        //also, we compute phi and psi here, node by node
        if (C != 0.0)
          {
          phi[nd0] = C*(alpha11*(w1x[nd0]/w1[nd0]) - beta*(w1y[nd0]/w1[nd0]))/alpha11;
          psi[nd0] = C*(-beta*(w2x[nd0]/w2[nd0]) + alpha12*(w2y[nd0]/w2[nd0]))/alpha12;
          }
        }

      //now, we place the appropriate weights in the appropriate entries of mat
      //notice that we will be filling in the row corresponding to the node in question and it's relationship to the other nodes in triangle, triangle by triangle
      //each node is related to the others in both directions, so we use += for the weights
      //notice, also, that forcing function have been included in off-diagonal nodes...multiplied by -2A since not divided by -2A!
      //w2* = w1*, but only for winslow; also, in winslow, the off diag elements are all zero, but this is taken care of by passing in coefficients
      //for L-E beta is set to zero
      w101 = alpha11*norm0a[0]*tx - 2.0*beta*norm0a[1]*tx + alpha12*norm0a[1]*ty;
      w101 /= -2.0*triarea;

      w102 = theta11*norm0a[1]*tx + theta12*norm0a[0]*ty;
      w102 /= -2.0*triarea;

      w111 = alpha11*norm1a[0]*tx - 2.0*beta*norm1a[1]*tx + alpha12*norm1a[1]*ty - 2.0*triarea*0.5*alpha11*phi[nd0]*tx;
      w111 /= -2.0*triarea;

      w112 = theta11*norm1a[1]*tx + theta12*norm1a[0]*ty;
      w112 /= -2.0*triarea;

      w121 = alpha11*norm2a[0]*tx - 2.0*beta*norm2a[1]*tx + alpha12*norm2a[1]*ty - 2.0*triarea*0.5*alpha12*psi[nd0]*ty;
      w121 /= -2.0*triarea;
   
      w122 = theta11*norm2a[1]*tx + theta12*norm2a[0]*ty;
      w122 /= -2.0*triarea;

      w201 = alpha21*norm0a[0]*tx - 2.0*beta*norm0a[1]*tx + alpha22*norm0a[1]*ty;
      w201 /= -2.0*triarea;

      w202 = theta21*norm0a[1]*tx + theta22*norm0a[0]*ty;
      w202 /= -2.0*triarea;

      w211 = alpha21*norm1a[0]*tx - 2.0*beta*norm1a[1]*tx + alpha22*norm1a[1]*ty - 2.0*triarea*0.5*alpha21*phi[nd0]*tx;
      w211 /= -2.0*triarea;

      w212 = theta21*norm1a[1]*tx + theta22*norm1a[0]*ty;
      w212 /= -2.0*triarea;

      w221 = alpha21*norm2a[0]*tx - 2.0*beta*norm2a[1]*tx + alpha22*norm2a[1]*ty - 2.0*triarea*0.5*alpha22*psi[nd0]*ty;
      w221 /= -2.0*triarea;
   
      w222 = theta21*norm2a[1]*tx + theta22*norm2a[0]*ty;
      w222 /= -2.0*triarea;

      /*//debug
      if (fabs(w101) < 1.0e-15 || fabs(w111) < 1.0e-15 || fabs(w121) < 1.0e-15 || fabs(w201) < 1.0e-15 || fabs(w211) < 1.0e-15 || fabs(w221) < 1.0e-15)
        printf("\nYou have a 0 weight on a submatrix diagonal at triangle %d",i);*/

      //now, use these to fill in mat!
      //fill in diagonal entries with w101, w201 since these u, v respectively
      //fill in off-diagonal entries with w102, w202 since these are for v, u respectively (causes coupling)
      //notice, in winslow, the off diag entries of the 2x2 remain 0
      mat[iau[nd0]][0][0] += w101; //u
      mat[iau[nd0]][0][1] += w102; //v
      mat[iau[nd0]][1][0] += w202; //u
      mat[iau[nd0]][1][1] += w201; //v
        
      //now, fill in off diag with n1
      //loop through n0's row to find column for n1
      for (k = ia[nd0]; k < ia[nd0+1]; k++)
        {
        if (ja[k] == nd1)
          {
          mat[k][0][0] += w111; //u
          mat[k][0][1] += w112; //v
          mat[k][1][0] += w212; //u
          mat[k][1][1] += w211; //v
          }
        }

      //now, fill in off diag with n2
      //loop through n0's row to find column for n2
      //notice, we can "go out of bounds" in ia since we added that extra element of the "final" node number (computed as the number of nodes)
      for (k = ia[nd0]; k < ia[nd0+1]; k++)
        {
        if (ja[k] == nd2)
          {
          mat[k][0][0] += w121; //u
          mat[k][0][1] += w122; //v
          mat[k][1][0] += w222; //u
          mat[k][1][1] += w221; //v
          }
        }
      } //end each node of triangle loop
    } //end triangle loop

  //debug
  //set high/low to start process
  double maxphi = -1.0e20;
  double maxpsi = -1.0e20;
  double minphi = 1.0e20;
  double minpsi = 1.0e20;

  for (i = 0; i < nn; i++)
    {
    maxphi = MAX(maxphi,phi[i]);
    minphi = MIN(minphi,phi[i]);
    maxpsi = MAX(maxpsi,psi[i]);
    minpsi = MIN(minpsi,psi[i]);
    }

  printf("\nMax phi = %lf\n",maxphi);
  printf("\nMax psi = %lf\n",maxpsi);
  printf("\nMin phi = %lf\n",minphi);
  printf("\nMin psi = %lf\n",minpsi);

  /*//debug
  flag = 0;
  for (i = 0; i < nn; i++)
    {
    if (fabs(mat[iau[i]][0][0]) < 1.0e-15 || fabs(mat[iau[i]][1][1]) < 1.0e-15)
      {
      printf("\n mat [%d] has a weight of 0 and the node is tagged %d!\n",i,tag[i]);
      flag++;
      }
    }
  printf("\n total 2x2s filled with zeros = %d!\n",flag);*/

  //now, we invert the diagonal matrices and place them right back in mat since we only need inverted versions to do solve
  for (i = 0; i < nn; i++)
    {
    //even though there are would-be placeholder variables (0,1) and (1,0) that are zeros, we want to be generic for L-E
    //plus, this saves us memory access time
    a = mat[iau[i]][0][0];
    b = mat[iau[i]][0][1];
    c = mat[iau[i]][1][0];
    d = mat[iau[i]][1][1];

    det = a*d - b*c; //take determinant
    //if (fabs(det) < 1.0e-15)
      //printf("\n determinant for mat [%d] = 0!\n",iau[i]);
    det = 1.0/det; //invert to multiply by entries in inverse

    //now, recreate diagonal entries as inverses
    mat[iau[i]][0][0] = det*d;
    mat[iau[i]][0][1] = -det*b;
    mat[iau[i]][1][0] = -det*c;
    mat[iau[i]][1][1] = det*a;
    }
   
  return;
  } 
 
