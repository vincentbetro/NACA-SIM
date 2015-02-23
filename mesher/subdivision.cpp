#include "subdivision.h"

//takes variable info, number of vars, number of functions to test on, which functions to test on, edge or grad based refinement toggle, tri connectivity, physical nodes, std deviation inputs for coarsening or refinement, node to node hash, nabor conn, boundry tagging, coarsening map, power to raise side length to for refinement function

//returns new node total

int subdivision(double **q, int nv, int test_funct, int **tri, int nt, double **nodep, int nn, int funct_type, double Cr[3], double Cc[3], double power[3], List **hash, int **nbr, int *tag, int *map, int ft)
  {
  int i, j, k, n, l; //counters
  double x, y; //node coord placeholders
  double favg = 0.0, stddev = 0.0, fr = 0.0, fc = 0.0, f = 0.0; //avg of refinement funct, standard dev of refinement funct, refinement threshold, coarsening threshold, refinement function itself
  int n0, n1, n2, nd0, nd1; //node placeholders for cleanliness
  Point p0, p1, p2; //points for clean vector creation
  Vector s0, s1, s2, r0, r1, gradphi0, gradphi1; //vectors for side length, normals for gradient, gradient vectors
  double len, len0, len1, len2; //lengths of side vectors
  int newnodes = nn; //counter for new node numbers, set to start at nn (cannot do in loop, else resets)
  double gamma = 1.4;

  //declare funct (1-D, since one at a time)
  double *funct = (double*)calloc(nn,sizeof(double));
 
  //declare grad arrays, cvareaG, although only used if funct_type == 0
  double *dfdy = NULL, *dfdx = NULL;
  double *cvareaG = NULL;

  //only allocate for dfdy, dfdx, cvareaG if funct_type == 0  
  if (funct_type == 0)
    {  
    dfdy = (double*)calloc(nn,sizeof(double));
    dfdx = (double*)calloc(nn,sizeof(double));
    cvareaG = (double*)calloc(nn,sizeof(double));
    //first, we compute areas of CV's for gradient making, since only needs done once..only if funct_type == 0
    //pass triangle, node, boundary info into area_calc...cvareaG and tag will come out filled in
    area_calc(cvareaG, tag, tri, nt, nodep, nn);
    }

  //now, we must compute each function at each node and store
  //in the future, we will use q at each node to compute diff functions (not necessarily nv of them, just more than one)
  //we do these one at a time
  for (k = 0; k < test_funct; k++)
    { 
    //now, we compute the test function at each node  
    for (i = 0; i < nn; i++)
      {
      //if we have nv == 0, then we are using analytic
      if (nv == 0)
        {  
        x = nodep[i][0];
        y = nodep[i][1];
    
        //store outcome of analytic in funct (in future, we can store outcomes of diff funct in arbitrarily numbered sections of funct)
        if (ft == 0)
          funct[i] = analytic_circ(x,y);
        if (ft == 1)
          funct[i] = analytic_curve(x,y);
        }
      else if ((nv == 4 && ft == 2) || (nv == 4 && ft == 5 && k == 0))
        {
        funct[i] = (gamma - 1.0)*(q[i][3] - 0.5*((q[i][1]*q[i][1] + q[i][2]*q[i][2])/q[i][0]));
        }
      else if ((nv == 4 && ft == 3) || (nv == 4 && ft == 5 && k == 1))
        {
        funct[i] = sqrt((q[i][1]*q[i][1] + q[i][2]*q[i][2])/(q[i][0]*q[i][0]));
        }
      else if ((nv == 4 && ft == 4) || (nv == 4 && ft == 5 && k == 2))
        {
        double p = (gamma - 1.0)*(q[i][3] - 0.5*((q[i][1]*q[i][1] + q[i][2]*q[i][2])/q[i][0]));
        double v = sqrt((q[i][1]*q[i][1] + q[i][2]*q[i][2])/(q[i][0]*q[i][0]));
        funct[i] = v / sqrt((gamma*p)/q[i][0]);
        }
      }
  
    //now, we compute the refinement function one of two ways, but doing the same basic operations
    //only define gradients if using that type of function
    if (funct_type == 0)
      {
      //we compute grads of all functions
      //now, we pass in cvareaG, tag of bd nodes, tri info, node info, the funct array, as well as dfdx, dfdy empty to come out full
      gradients(cvareaG, tag, tri, nt, nodep, nn, funct, dfdx, dfdy);
      }

    favg = 0.0; //be sure to reset, just in case, since summing (even though done in declaration)
    stddev = 0.0; //be sure to reset, just in case, since summing (even though done in declaration)
    n = 0; //reset n to be sure we are summing edges from 0
  
    //now, we loop use node_node hash and determine favg (this means each edge (including bd) will be seen twice, making avg equal)
    for (i = 0; i < nn; i++)
      {
      for (j = 0; j < hash[i]->max; j++)
        {
        if (hash[i]->list[j] == i)
          continue;

        //now, we make edges out of all other nodes in hash

        //use node indices for cleanliness
        n0 = i;
        n1 = hash[i]->list[j];

        //use Points for clean vectors (unit and length mag)
        p0 = Point(nodep[n0][0],nodep[n0][1]);
        p1 = Point(nodep[n1][0],nodep[n1][1]);

        //create vectors along side
        s0 = Vector(p0,p1);

        //find its length
        len = s0.magnitude();

        //normalize it
        s0.normalize();
 
        //determine grad vectors and normals on edge
        if (funct_type == 0)
          {
          gradphi0 = Vector(dfdx[n0],dfdy[n0]);
          gradphi1 = Vector(dfdx[n1],dfdy[n1]);
          r0 = s0;
          r1 = Vector(-s0[0],-s0[1]);
          }           
          
        //now, compute f at both nodes, take average, summing whole time toward favg (and summing n to get appropriate factor to divide by)
        //be sure to compute based on user input
        if (funct_type == 0)
          favg += (fabs((gradphi0*r0))*pow(len,power[k]) + fabs((gradphi1*r1))*pow(len,power[k]))/2.0;
        if (funct_type == 1)
          favg += fabs((funct[n1] - funct[n0])/len)*pow(len,power[k]);
        n++; //sum number of edges
        }
      }

    //find favg now
    favg /= n;
    printf("\nfavg = %16.10e\n",favg);

    //now, we loop use node_node hash and determine stddev (this means each edge (including bd) will be seen twice, making deviation balanced)
    for (i = 0; i < nn; i++)
      {
      for (j = 0; j < hash[i]->max; j++)
        {
        if (hash[i]->list[j] == i)
          continue;

        //now, we make edges out of all other nodes in hash

        //use node indices for cleanliness
        n0 = i;
        n1 = hash[i]->list[j];

        //use Points for clean vectors (unit and length mag)
        p0 = Point(nodep[n0][0],nodep[n0][1]);
        p1 = Point(nodep[n1][0],nodep[n1][1]);

        //create vectors along side
        s0 = Vector(p0,p1);

        //find its length
        len = s0.magnitude();

        //normalize it
        s0.normalize();

        //determine grad vectors and normals on edge
        if (funct_type == 0)
          {
          gradphi0 = Vector(dfdx[n0],dfdy[n0]);
          gradphi1 = Vector(dfdx[n1],dfdy[n1]);
          r0 = s0;
          r1 = Vector(-s0[0],-s0[1]);
          }           
          
        //now, compute f at both nodes, take std dev, summing whole time toward stddev (n already summed from previous loop)
        //be sure to compute based on user input
        if (funct_type == 0)
          stddev += (((fabs((gradphi0*r0))*pow(len,power[k]) + fabs((gradphi1*r1))*pow(len,power[k]))/2.0) - favg)*(((fabs((gradphi0*r0))*pow(len,power[k]) + fabs((gradphi1*r1))*pow(len,power[k]))/2.0) - favg);
        if (funct_type == 1)
          stddev += ((fabs((funct[n1] - funct[n0])/len)*pow(len,power[k])) - favg)*((fabs((funct[n1] - funct[n0])/len)*pow(len,power[k])) - favg);
        }
      }

    //find final stddev now
    stddev /= n;
    stddev = sqrt(stddev);

    printf("\nStddev = %16.10e\n",stddev);

    //finally, determine fr and fc
    fr = favg + Cr[k]*stddev;
    printf("\nfr = %16.10e\n",fr);
    fc = favg - Cc[k]*stddev;
    printf("\nfc = %16.10e\n",fc);     

    //now, we loop through all triangles and determine whether or not to refine/delete an edge
    //we will also go ahead and create new nodes and inform nabors simultaneously
    //newnodes = nn; //set starting point for new node creation to current num of nodes (done externally, now)
    for (i = 0; i < nt; i++)
      {
      //use node indices for cleanliness
      n0 = tri[i][0];
      n1 = tri[i][1];
      n2 = tri[i][2];

      //use Points for clean vectors (unit and length mag)
      p0 = Point(nodep[n0][0],nodep[n0][1]);
      p1 = Point(nodep[n1][0],nodep[n1][1]);
      p2 = Point(nodep[n2][0],nodep[n2][1]);

      //create vectors along each side
      s0 = Vector(p0,p1);
      s1 = Vector(p1,p2);
      s2 = Vector(p2,p0);

      //find their length
      len0 = s0.magnitude();
      len1 = s1.magnitude(); 
      len2 = s2.magnitude();

      //normalize them
      s0.normalize();
      s1.normalize();
      s2.normalize();

      for (j = 0; j < 3; j++)
        {
        switch(j)
          {
          case 0:
            if (funct_type == 0)
              {
              gradphi0 = Vector(dfdx[n0],dfdy[n0]);
              gradphi1 = Vector(dfdx[n1],dfdy[n1]);
              r0 = s0;
              r1 = Vector(-s0[0],-s0[1]); 
              }
            len = len0;          
            nd0 = n0;
            nd1 = n1;  
          break;
          case 1:
            if (funct_type == 0)
              {
              gradphi0 = Vector(dfdx[n1],dfdy[n1]);
              gradphi1 = Vector(dfdx[n2],dfdy[n2]);
              r0 = s1;
              r1 = Vector(-s1[0],-s1[1]); 
              }
            len = len1;          
            nd0 = n1;
            nd1 = n2;
          break;
          case 2:
            if (funct_type == 0)
              {
              gradphi0 = Vector(dfdx[n2],dfdy[n2]);
              gradphi1 = Vector(dfdx[n0],dfdy[n0]);
              r0 = s2;
              r1 = Vector(-s2[0],-s2[1]); 
              }
            len = len2;          
            nd0 = n2;
            nd1 = n0;          
          break;
          default:
            printf("\nYou have a triangle with more than three sides according to the loop counter in subdivison.cpp!\n");
            fflush(stdout);
            exit(0);
          break;
          }
          
        //now, find f again, and determine if it is within the refinement threshold
        //be sure to compute based on user input
        if (funct_type == 0)
          f = (fabs((gradphi0*r0))*pow(len,power[k]) + fabs((gradphi1*r1))*pow(len,power[k]))/2.0;
        if (funct_type == 1)
          f = fabs((funct[nd1] - funct[nd0])/len)*pow(len,power[k]);

        //since we only need one function to represent for refinement, we create new node (if already created, forget about it!)
        if (f  > fr && tri[i][j+3] < 0)
          {
          tri[i][j+3] = newnodes;
          //also, we need to notify nabor (if one exists)
          n = nbr[i][j]; //locate nabor
          //if we are on a boundary, this is not an issue, and we skip adding node to nabor
          if (n >= 0)
            {
            //locate node it shares with tri[i]
            for (l = 0; l < 3; l++)
              {
              if (nbr[n][l] == i)
                tri[n][l+3] = newnodes; //reset its node number
              }
            }           
          newnodes++; //increment newnodes
          }

        //notice that we will mark edge nodes for deletion thru map, since ALL functions must agree
        //we use reverse logic, so if we are greater than fc, we mark to keep, thus if it passes thru all functions without being kept, we delete it
	//the boundary nodes have all been reset to zero, so we will automatically NOT delete them!
        if (f > fc)
          {
          if (map[nd0] == -1)
            map[nd0] = 0;
          if (map[nd1] == -1)
            map[nd1] = 0;
          }
	//we do not want to change newnodes here...we wait to delete until after all reallocation and renumbering done
        }
      }
    }

  //clean up memory
  freenull(funct);
  freenull(dfdy);
  freenull(dfdx);
  freenull(cvareaG);

  return(newnodes);
  } 
