#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Linked_List.h"
#include "List.h"
#include "Lhash.h"
#include "crs.h"
#include "winslow.h"
#include "linear_elastic.h"
#include "trimesh.h"
#include "Point.h"
#include "Vector.h"
#include "make_nbrs_HO.h"
#include "subdivision.h"
#include "refine_tri.h"
#include "opt_smooth.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))
#define freenull(a) {if(a) free(a); a=NULL;}

//Generally, run to cvg of 1.0e-6

int main(int argcs, char* pArgs[])
{
  int nn = 0, nblk = 0, nt = 0, nq = 0, nb = 0, nv = 0, nc = 0; //number of nodes, number of blocks, number of triangles, number of quads, number of boundaries, num vars, num const
  int nnp, nblkp, ntp, nqp, nbp; //p signifies physical mesh version to be checked against above
  //placeholders to read physical mesh info into to check against computational (saves storing in more matrices and wasting memory)
  //also, flag for derefining (set to zero for freeing mem) and for remeshing (set to zero for freeing mem) and for reallocing tri within trimesh loop
  int trip0, trip1, trip2, quadp0, quadp1, quadp2, quadp3, bsp0, bsp1, nbsp, flagC = 0, remesh = 0, tdim; 
  int i, j, k, n, m, l, s; //counters
  int old_nn; //old number of nodes, needed for freeing memory if operation == 1
  const int bdim = 132; //buffer size
  char buff[bdim]; //buffer
  FILE *fp = NULL; //file pointer to be used for both files
  //declare matrices for x,y coords of computational/physical nodes (indexed), triangle conn, quad conn
  double **nodec = NULL, **nodep = NULL;
  int **tri = NULL, **quad = NULL; 
  //boundary segments, indexed as boundary number, edge number, node number on edge; number boundary segments per edge 
  int ***bs = NULL, *nbs = NULL;
  //array for each node's values of variables, array for any constants
  double **q = NULL, *constants = NULL; //NOTE:  These are read from physical mesh ONLY   
  
  printf("\n====================================================");
  printf("\n   Driver for 2-D Winslow/Linear-Elastic Smoother   ");
  printf("\n====================================================\n");

  //be sure both physical and computational are supplied, since connectivity must be identical, while node placement needn't be
  if (--argcs < 2)
    {
    printf("\nYou must specify both a computational mesh and physical mesh file, in that order!");
    fflush(stdout);
    exit(0);
    }

  //check to be sure computational mesh file can be opened
  if ((fp=fopen(pArgs[argcs],"r")) == 0)
    {
    printf("\nCouldn't open computational mesh file <%s>\n",pArgs[argcs]);
    fflush(stdout);
    exit(0);
    }

  // read number of nodes
  fgets(buff,bdim,fp); //skip commment
  fgets(buff,bdim,fp); //read nn
  sscanf(buff,"%d",&nn);
  printf("\nNumber of nodes = %d",nn);
  fflush(stdout); //be sure prints to screen

  // allocate for computational nodes, using calloc, due to possible realloc
  nodec = (double**)calloc(nn,sizeof(double*));
  for (i = 0; i < nn; i++) 
    nodec[i] = (double*)calloc(2,sizeof(double));

  // read in coordinates
  for (i = 0; i < nn; i++)
    {
    fgets(buff,bdim,fp); //read node coords to be placed in row of node number (c-indexed)
    sscanf(buff,"%lf %lf",&(nodec[i][0]),&(nodec[i][1]));
    }

  // read number of blocks
  fgets(buff,bdim,fp); //skip commment
  fgets(buff,bdim,fp); //read nblk
  sscanf(buff,"%d",&nblk);
  printf("\nNumber of blocks = %d",nblk);
  fflush(stdout); //be sure prints to screen

  // read number of triangles
  fgets(buff,bdim,fp); //skip commment
  fgets(buff,bdim,fp); //read nt
  sscanf(buff,"%d",&nt);
  printf("\nNumber of triangles = %d",nt);
  fflush(stdout); //be sure prints to screen

  // allocate for triangle connectivity, using calloc, due to possible realloc
  tri = (int**)calloc(nt,sizeof(int*));
  for (i = 0; i < nt; i++) 
    tri[i] = (int*)calloc(3,sizeof(int));

  // read in triangle connectivity
  for (i = 0; i < nt; i++)
    {
    fgets(buff,bdim,fp); //read ccw nodes for each triangle
    sscanf(buff,"%d %d %d",&(tri[i][0]),&(tri[i][1]),&(tri[i][2]));
    }
  
  // decrement tri by 1 to make fortran-indexed into c-indexed
  for (i = 0; i < nt; i++)
    for (j = 0; j < 3; j++)
      tri[i][j]--;

  // read number of quads
  fgets(buff,bdim,fp); //skip commment
  fgets(buff,bdim,fp); //read nq
  sscanf(buff,"%d",&nq);
  printf("\nNumber of quads = %d",nq);
  fflush(stdout); //be sure prints to screen

  if (nq > 0)
    {
    // allocate for quad connectivity, using calloc, due to possible realloc
    quad = (int**)calloc(nq,sizeof(int*));
    for (i = 0; i < nq; i++) 
      quad[i] = (int*)calloc(4,sizeof(int));

    // read in quad connectivity
    for (i = 0; i < nq; i++)
      {
      fgets(buff,bdim,fp); //read ccw nodes for each quad
      sscanf(buff,"%d %d %d %d", &(quad[i][0]), &(quad[i][1]), &(quad[i][2]), &(quad[i][3]));
      }

    // decrement quad by 1 to make fortran-indexed into c-indexed
    for (i = 0; i < nq; i++)
      for (j = 0; j < 4; j++)
        quad[i][j]--;
    }

  // read number of boundaries
  fgets(buff,bdim,fp); //skip commment
  fgets(buff,bdim,fp); //read nb
  sscanf(buff,"%d",&nb);
  printf("\nNumber of boundaries = %d",nb);
  fflush(stdout); //be sure prints to screen

  // allocate for boundary segment list (first index only...boundary number) and number of segments per boundary, using calloc, due to possible realloc
  // 3-D array so contiguous not needed
  bs = (int***)calloc(nb,sizeof(int**));
  nbs = (int*)calloc(nb,sizeof(int));

  // read number of edges per boundary as well as corresponding node numbers in a loop
  for (i = 0; i < nb; i ++)
    {
    //read in number of edges
    fgets(buff,bdim,fp); //skip commment
    fgets(buff,bdim,fp); //read nbs[i] (i is boundary number)
    sscanf(buff,"%d",&(nbs[i]));
    printf("\nNumber of edges on boundary %d = %d",i,nbs[i]);
    fflush(stdout); //be sure prints to screen

    //now, allocate this many edges in the area of bs corresponding to boundary i
    bs[i] = (int**)calloc(nbs[i],sizeof(int*));

    //now, allocate two spots per edge (for node numbers) in final dimension of bs
    for (j = 0; j < nbs[i]; j++)
      bs[i][j] = (int*)calloc(2,sizeof(int));
   
    //now, read in node numbers
    for (j = 0; j < nbs[i]; j++)
      {
      fgets(buff,bdim,fp); //read nodes of each edge
      sscanf(buff,"%d %d", &(bs[i][j][0]), &(bs[i][j][1]));
      }
    
    // now, decrement edges by 1 to make fortran-indexed into c-indexed
    for (j = 0; j < nbs[i]; j++)
      for (k = 0; k < 2; k++)
        bs[i][j][k]--;
    }

  //we do not need further data from file, so we close and prepare to read physical mesh and compare its connectivity to computational mesh
  fclose(fp);

  //check to be sure physical mesh file can be opened
  if ((fp=fopen(pArgs[--argcs],"r")) == 0)
    {
    printf("\nCouldn't open physical mesh file <%s>\n",pArgs[argcs]);
    exit(0);
    }

  // read number of nodes
  fgets(buff,bdim,fp); //skip commment
  fgets(buff,bdim,fp); //read nnp
  sscanf(buff,"%d",&nnp);
  printf("\nNumber of physical mesh nodes = %d",nnp);
  fflush(stdout); //be sure prints to screen
  //check for agreement between physical and computational meshes
  if (nnp != nn)
    {
    printf("\nYour physical and computational mesh number of nodes do not agree!");
    fflush(stdout);
    exit(0);
    }

  // allocate for physical nodes, using calloc, due to possible realloc
  nodep = (double**)calloc(nnp,sizeof(double*)); 
  for (i = 0; i < nnp; i++) 
    nodep[i] = (double*)calloc(2,sizeof(double));

  // read in coordinates
  for (i = 0; i < nnp; i++)
    {
    fgets(buff,bdim,fp); //read node coords to be placed in row of node number (c-indexed)
    sscanf(buff,"%lf %lf",&(nodep[i][0]),&(nodep[i][1]));
    }

  // read number of blocks
  fgets(buff,bdim,fp); //skip commment
  fgets(buff,bdim,fp); //read nblkp
  sscanf(buff,"%d",&nblkp);
  printf("\nNumber of physical mesh blocks = %d",nblkp);
  fflush(stdout); //be sure prints to screen
  //check for agreement between physical and computational meshes
  if (nblkp != nblk)
    {
    printf("\nYour physical and computational mesh number of blocks do not agree!");
    fflush(stdout);
    exit(0);
    }

  // read number of triangles
  fgets(buff,bdim,fp); //skip commment
  fgets(buff,bdim,fp); //read ntp
  sscanf(buff,"%d",&ntp);
  printf("\nNumber of physical mesh triangles = %d",ntp);
  fflush(stdout); //be sure prints to screen
  //check for agreement between physical and computational meshes
  if (ntp != nt)
    {
    printf("\nYour physical and computational mesh number of triangles do not agree!");
    fflush(stdout);
    exit(0);
    }

  // read in triangle connectivity
  for (i = 0; i < ntp; i++)
    {
    fgets(buff,bdim,fp); //read ccw nodes for each triangle
    sscanf(buff,"%d %d %d",&trip0,&trip1,&trip2);
    //decrement trip0, trip1, trip2 for comparison in if statement
    if (trip0 - 1 != tri[i][0] || trip1 - 1 != tri[i][1] || trip2 - 1 != tri[i][2])
      {
      printf("\nYour physical and computational mesh triangle connectivities do not agree!");
      fflush(stdout);
      exit(0);
      }
    }

  // read number of quads
  fgets(buff,bdim,fp); //skip commment
  fgets(buff,bdim,fp); //read nqp
  sscanf(buff,"%d",&nqp);
  printf("\nNumber of physical mesh quads = %d",nqp);
  fflush(stdout); //be sure prints to screen
  //check for agreement between physical and computational meshes
  if (nqp != nq)
    {
    printf("\nYour physical and computational mesh number of quads do not agree!");
    fflush(stdout);
    exit(0);
    }

  if (nqp > 0)
    {
    // read in quad connectivity
    for (i = 0; i < nqp; i++)
      {
      fgets(buff,bdim,fp); //read ccw nodes for each quad, decrement by 1 to make fortran-indexed, c-indexed
      sscanf(buff,"%d %d %d %d",&quadp0,&quadp1,&quadp2,&quadp3);
      //decrement quadp0, quadp1, quadp2, quadp3 for comparison in if statement
      if (quadp0 - 1 != quad[i][0] || quadp1 - 1 != quad[i][1] || quadp2 - 1 != quad[i][2] || quadp3 - 1 != quad[i][3])
        {
        printf("\nYour physical and computational mesh quad connectivities do not agree!");
        fflush(stdout);
        exit(0);
        }
      }
    }

  // read number of boundaries
  fgets(buff,bdim,fp); //skip commment
  fgets(buff,bdim,fp); //read nbp
  sscanf(buff,"%d",&nbp);
  printf("\nNumber of physical mesh boundaries = %d",nbp);
  fflush(stdout); //be sure prints to screen
  //check for agreement between physical and computational meshes
  if (nbp != nb)
    {
    printf("\nYour physical and computational mesh number of boundaries do not agree!");
    fflush(stdout);
    exit(0);
    }

  // read number of edges per boundary as well as corresponding node numbers in a loop
  for (i = 0; i < nbp; i ++)
    {
    //read in number of edges
    fgets(buff,bdim,fp); //skip commment
    fgets(buff,bdim,fp); //read nbsp
    sscanf(buff,"%d",&nbsp);
    printf("\nNumber of edges on boundary %d in physical mesh = %d",i,nbsp);
    fflush(stdout); //be sure prints to screen
    //check for agreement between physical and computational meshes
    if (nbsp != nbs[i])
      {
      printf("\nYour physical and computational mesh number of edges per boundary do not agree!");
      fflush(stdout);
      exit(0);
      }
   
    //now, read in node numbers, again switching form fortran to c indexing
    for (j = 0; j < nbsp; j++)
      {
      fgets(buff,bdim,fp); //read nodes of each edge
      sscanf(buff,"%d %d", &bsp0, &bsp1);
      //decrement bs0, bs1 for comparison in if statement
      if (bsp0 - 1 != bs[i][j][0] || bsp1 - 1 != bs[i][j][1])
        {
        printf("\nYour physical and computational mesh boundary edge indices do not agree!");
        fflush(stdout);
        exit(0);
        }
      }
    }

  //NOTE:  If using old mesh files without variable info, comment out!
  // read number of constants
  fgets(buff,bdim,fp); //skip commment
  fgets(buff,bdim,fp); //read nc
  sscanf(buff,"%d",&nc);
  printf("\nNumber of constants = %d",nc);
  fflush(stdout); //be sure prints to screen

  if (nc > 0)
    {
    // allocate for constants
    constants = (double*)calloc(nc,sizeof(double)); //allocate
    
    // read in constants (we will need to know what each is in order, a priori)
    for (i = 0; i < nc; i++)
      {
      fgets(buff,bdim,fp); //read constant
      sscanf(buff,"%lf", &(constants[i]));
      }
    }

  // read number of vars (usually four, but we could code for more if needed)
  fgets(buff,bdim,fp); //skip commment
  fgets(buff,bdim,fp); //read nv
  sscanf(buff,"%d",&nv);
  printf("\nNumber of variables per node = %d",nv);
  fflush(stdout); //be sure prints to screen
  for (i = 0; i < nv; i++)
    fgets(buff,bdim,fp); //skip nv commments on what they are named

  if (nv > 0 && nv < 5)
    {
    // allocate for each nodes q values
    q = (double**)calloc(nn,sizeof(double*));
    for (i = 0; i < nn; i++) 
      q[i] = (double*)calloc(nv,sizeof(double));

    // read in q at each node
    for (i = 0; i < nn; i++)
      {
      fgets(buff,bdim,fp); //read q in order
      sscanf(buff,"%lf %lf %lf %lf", &(q[i][0]), &(q[i][1]), &(q[i][2]), &(q[i][3]));
      }
    }

  if (nv > 4)
    {
    printf("\nThis code is 2-D and thus only set up for density, x-momentum, y-momentum, energy.  Exiting....\n");
    fflush(stdout);
    exit(0);
    }

  //we do not need further data from file, so we close
  fclose(fp);

  //now, read in desired output file name
  char filename[bdim]; //to hold user input
  printf("\nEnter desired name of physical mesh, gnuplot file, without file extension ->");
  fgets(buff,bdim,stdin); //read in filename
  sscanf(buff,"%s",&filename); //place in variable
  int n0, n1, n2; //node indices to keep things clean and save memory access
  
  //read in operation to be done, smoothing/movement or subdivision refinement
  int operation = 0; //init to smoothing/mesh movement
  printf("\nEnter operation to be done: 0 for smoothing/mesh movement OR 1 for subdivision refinement  ->");
  fgets(buff,bdim,stdin); //read in operation
  sscanf(buff,"%d",&operation); //place in variable
  //error checking
  if (operation < 0 || operation > 1)
    {
    printf("\nOperation must be 0 for smoothing/mesh movement or 1 for subdivision refinement.  Exiting....");
    fflush(stdout);
    exit(0);
    }
  //print out op to be done
  if (operation == 0)
    {
    printf("\nPerforming smoothing/mesh movement. Continuing....");
    fflush(stdout);
    }
  if (operation == 1)
    {
    printf("\nPerforming subdivision refinement. Continuing....");
    fflush(stdout);
    }

  if (operation == 0)
    {
    //prepare to read in information about omega, inner iter, outer iter, and convergence level
    int iterin = 0, iterout = 0; //number of iterations to perform on solve before recomputing alpha, gamma, beta, gradients/ cutoff iterations to perform if cvg not reached sooner
    double omega, cvg; //relaxation factor/ desired level of convergence
    int smooth_type; //0 for winslow, 1 for L-E, 2 for optimization based
    double poisson = 0.0; //poisson's ratio, init to zero since will only be set if smooth_type = 1
    double time, time_slice, alphamax; //read in values used to compute delta t and angle airfoil to be moved through
    int bd1, bd2; //holds airfoil boundry numbers (to be moved), C++ indexed
    double A, B, C; //forcing function coefficients
    int whichff; //decide which forcing function to use

    //read in type of smoothing to be done, if L-E, read in poisson's ratio
    printf("\nEnter desired type of smoothing: 0 for Winslow, 1 for Linear-Elastic, 2 for optimization-based  ->");
    fgets(buff,bdim,stdin); //read in smooth_type
    sscanf(buff,"%d",&smooth_type); //place in variable
    //error checking
    if (smooth_type == 0)
      {
      printf("\nWill perform Winslow Smoothing.  Continuing....");
      //also, read in inner and outer iterative loop limits here, since L-E only has inner loop
      printf("\nEnter number of iterations for inner point-iterative GS-SSOR loop ->");
      fgets(buff,bdim,stdin); //read in iterin
      sscanf(buff,"%d",&iterin); //place in variable
      //error checking
      if (iterin <= 0)
	{
	printf("\nCondition to satisfy: iterin > 0.  Exiting....");
        fflush(stdout);
	exit(0);
	}
      printf("\nEnter number of iterations for outer point-iterative GS-SSOR loop ->");
      fgets(buff,bdim,stdin); //read in iterout
      sscanf(buff,"%d",&iterout); //place in variable
      //error checking
      if (iterout <= 0)
	{
	printf("\nCondition to satisfy: iterout > 0.  Exiting....");
        fflush(stdout);
	exit(0);
	}
      //get forcing function constants
      printf("\nEnter first degree constant for forcing functions (A) ->");
      fgets(buff,bdim,stdin); //read in A
      sscanf(buff,"%lf",&A); //place in variable
      printf("\nEnter second degree constant for forcing functions (B) ->");
      fgets(buff,bdim,stdin); //read in B
      sscanf(buff,"%lf",&B); //place in variable
      printf("\nEnter general constant for forcing functions (C) (for none, use 0.0) ->");
      fgets(buff,bdim,stdin); //read in C
      sscanf(buff,"%lf",&C); //place in variable

      //get type of forcing function
      printf("\nEnter 0 to use pressure, 1 for velocity magnitude, or 2 for mach number ->");
      fgets(buff,bdim,stdin); //read in whichff
      sscanf(buff,"%d",&whichff); //place in variable
      //error checking
      if (whichff < 0 || whichff > 2)
	{
	printf("\nYou must enter 0 to use pressure, 1 for velocity magnitude, or 2 for mach number.  Exiting....");
        fflush(stdout);
	exit(0);
	}
      }
    else if (smooth_type == 1)
      {
      printf("\nWill perform Linear-Elastic Smoothing.  Continuing....");
      printf("\nEnter Poisson's Ratio (0.0 < nu < 0.5)  ->");
      fgets(buff,bdim,stdin); //read in poisson
      sscanf(buff,"%lf",&poisson); //place in variable
      //error checking
      if (poisson <= 0.0 || poisson >= 0.5)
	{
	printf("\nPoisson's ratio MUST be between 0.0 and 0.5!  Exiting....");
        fflush(stdout);
	exit(0);
	}
      //read in iterative slove loop limit here (only "inner" loop for L-E)
      printf("\nEnter number of iterations for point-iterative GS-SSOR loop ->");
      fgets(buff,bdim,stdin); //read in iterin
      sscanf(buff,"%d",&iterin); //place in variable
      //error checking
      if (iterin <= 0)
	{
	printf("\nCondition to satisfy: iterin > 0.  Exiting....");
        fflush(stdout);
	exit(0);
	}
      }
    else if (smooth_type == 2)
      {
      printf("\nWill perform Optimization-Based Smoothing.  Continuing....");
      //read in convergence limit here
      printf("\nEnter number of iterations to be performed ->");
      fgets(buff,bdim,stdin); //read in iterin
      sscanf(buff,"%d",&iterin); //place in variable
      //error checking
      if (iterin <= 0)
	{
	printf("\nCondition to satisfy: iterations > 0.  Exiting....");
        fflush(stdout);
	exit(0);
	}
      }
    else
      {
      printf("\nYou must choose 0 for Winslow, 1 for Linear-Elastic, or 2 for Optimization-Based!  Exiting....");
      fflush(stdout);
      exit(0);
      }
    //read in info to do moving airfoil
    //read in amount of time to move airfoil
    printf("\nEnter amount of time to move the airfoil (seconds)  ->");
    fgets(buff,bdim,stdin); //read in time
    sscanf(buff,"%lf",&time); //place in variable
    //error checking
    if (time <= 0.0)
      {
      printf("\nPlease use a positive time in seconds greater than 0!  Exiting....");
      fflush(stdout);
      exit(0);
      }
    //read in number of time slices
    printf("\nEnter number of time slices at which to perform smoothing  ->");
    fgets(buff,bdim,stdin); //read in time_slice
    sscanf(buff,"%lf",&time_slice); //place in variable
    //error checking
    if (time_slice <= 0.0)
      {
      printf("\nPlease use a positive time slice greater than 0!  Exiting....");
      fflush(stdout);
      exit(0);
      }
    //read in degrees to be moved through over two cycles from 0 to alphamax
    printf("\nEnter degrees to rotate airfoil (per cycle for two cycles) over specified time period  ->");
    fgets(buff,bdim,stdin); //read in alphamax
    sscanf(buff,"%lf",&alphamax); //place in variable
    //error checking
    //if (alphamax < 0.0 || alphamax >= 360.0)
      //{
      //printf("\nPlease use a positive angle, between 0 (inclusive) and 360 degrees!  Exiting....");
      //fflush(stdout);
      //exit(0);
      //}
    //read in two boundaries we wish to move
    printf("\nEnter two boundaries that will be moved, C++ indexed, with a space in between  ->");
    fgets(buff,bdim,stdin); //read in bd1, bd2
    sscanf(buff,"%d %d",&bd1,&bd2); //place in variables
    //error checking
    if (bd1 >= nb)
      {
      printf("\nBoundary %d does not exist, since nb = %d.  Check your indexing!  Exiting....",bd1,nb);
      fflush(stdout);
      exit(0);
      }
    if (bd2 >= nb)
      {
      printf("\nBoundary %d does not exist, since nb = %d.  Check your indexing!  Exiting....",bd2,nb);
      fflush(stdout);
      exit(0);
      }
    //only need this info if doing winslow or L-E    
    if (smooth_type == 0 || smooth_type == 1)
      {
      //read in omega
      printf("\nEnter relaxation factor (omega) ->");
      fgets(buff,bdim,stdin); //read in omega
      sscanf(buff,"%lf",&omega); //place in variable
      //error checking
      if (omega <= 0.0 || omega > 1.0)
        {
        printf("\nCondition to satisfy: 0.0 < omega <= 1.0.  Exiting....");
        fflush(stdout);
        exit(0);
        }
      }
    //read in cvg value (even though doubtful to reach with Opt Based) 
    printf("\nEnter desired convergence level ->");
    fgets(buff,bdim,stdin); //read in cvg
    sscanf(buff,"%lf",&cvg); //place in variable
    //error checking
    if (cvg <= 1.0e-15 || cvg >= 1.0)
      {
      printf("\nCondition to satisfy: 1.0e-15 < cvg < 1.0.  Exiting....");
      fflush(stdout);
      exit(0);
      }

    //allocate for Lhash to be passed into routine and filled in
    List **hash;
    hash = new List*[nn]; //each node has a list of nodes it's attached to
    for (i=0; i < nn; i++)
      hash[i] = new List(); //make each node's list

    //now, create hash table to be passed to compressed row storage routine
    Lhash(nn, nt, tri, hash);

    //now, we build a node 2 cell connectivity
    List **NChash;

    //first, initialize node to element hash table
    NChash = new List*[nn];

    for (i = 0; i < nn; i++)
      NChash[i] = new List();

    //now, create node to element hash
    for (i = 0; i < nt; i++)
      for (j = 0; j < 3; j++)
        NChash[tri[i][j]]->Check_List(i);

    /*char *tester = "nodecellhash.dat";
    printf("\nOutput Hash File Filename = <%s>\n",tester);
    // Open file for write
    if ((fp = fopen(tester,"w")) == 0)
      {
      printf("\nError opening file <%s>.",tester);
      exit(0);
      }
    for (i=0; i < nn; i++)
      {
      for (j = 0; j < NChash[i]->max; j++)
        {
        fprintf(fp,"%d ",NChash[i]->list[j]);
        }
      fprintf(fp,"\n");
      }
    
    fclose(fp);*/


    //only do this for winslow, L-E, but declare vars external
    int mdim; //number of entries in hash table
    int *ia, *ja, *iau; //compressed row storage matrices
    double ***mat; //holds 2x2 weight matrices for each node
    if (smooth_type == 0 || smooth_type == 1)
      {
      //now, determine sizes of matrices to be passed into compressed row storage empty and passed out full
      mdim = 0; 
      for (i = 0; i < nn; i++)
        {
        mdim+=hash[i]->max;
        }

      //printf("\nmdim = %d",mdim);

      //now, allocate for all compressed row storage matrices     
      ia = (int*)calloc(nn+1,sizeof(int)); //ia will hold index (in ja) of starting point for each node (row) in smoothing matrix (mat)
      ja = (int*)calloc(mdim,sizeof(int)); //ja holds entries of hash table in row major order
      iau = (int*)calloc(nn,sizeof(int)); //iau holds index (in ja) of the diagonal elements in smoothing matrix (mat)

      //now, create compressed row storage version of matrix to be used to do smoothing
      crs(nn, mdim, ia, ja, iau, hash);

      //after creating compressed row storage, hash is no longer needed and can be deleted to save memory
      for (i=0; i < nn; i++)
        delete hash[i];
      delete[] hash;

      //now, allocate for mat (holds 2x2 weight matrices for each node)
      mat = (double***)calloc(mdim,sizeof(double**));
      for (i = 0; i < mdim; i++)
        {
        mat[i] = (double**)calloc(2,sizeof(double*));
        for (j = 0; j < 2; j++)
	  {
	  mat[i][j] = (double*)calloc(2,sizeof(double));
	  }
        }
      } //end only winslow, l-e

    //now, we need to tag boundary nodes so we don't do linear solve on them (in L-E we also use to create delu and delv)
    int *tag; 
    //create a tag array, interior = 1, boundary = 0
    tag = (int*)calloc(nn,sizeof(int));

    //initialize to one, and go back and make boundary nodes 0
    for (i = 0; i < nn; i++)
      tag[i] = 1;
    
    //now, loop through boundary nodes and tag them 0
    for (i = 0; i < nb; i++)
      for (j = 0; j < nbs[i]; j++)
	for (k = 0; k < 2; k++)
	    tag[bs[i][j][k]] = 0; //retag dups since still just a memory access and not doing so disallows cache coherence

    //we need to make gnuplot of original first
    //concatenate
    strcat(filename,"_00.dat");
    printf("\nOutput Plot File Filename = <%s>\n",filename);
    // Open file for write
    if ((fp = fopen(filename,"w")) == 0)
      {
      printf("\nError opening file <%s>.",filename);
      exit(0);
      }
    for (i=0; i < nt; i++)
      {
      n0 = tri[i][0];
      n1 = tri[i][1];
      n2 = tri[i][2];
      fprintf(fp,"%19.10e %19.10e 0.0\n",  nodep[n0][0],nodep[n0][1]);
      fprintf(fp,"%19.10e %19.10e 0.0\n",  nodep[n1][0],nodep[n1][1]);
      fprintf(fp,"%19.10e %19.10e 0.0\n",  nodep[n2][0],nodep[n2][1]);
      fprintf(fp,"%19.10e %19.10e 0.0\n\n",nodep[n0][0],nodep[n0][1]);
      }
    
      fclose(fp);

    //also, let's make a temp array of nodep, so we can find RMS error between orig physical and final physical (ostensibly in same location)
    //notice that winslow should recover while L-E will not recover exactly
    double **tnode;
    tnode = (double**)calloc(nn,sizeof(double*));
    for (i = 0; i < nn; i++)
      tnode[i] = (double*)calloc(2,sizeof(double));
    //init to orig nodep
    for (i = 0; i < nn; i++)
      {
      tnode[i][0] = nodep[i][0];
      tnode[i][1] = nodep[i][1];
      }

    //now, we prepare to rotate our airfoil with linear elastic or winslow smoothing
    //notice how we will use smooth_type to determine when and where we use physical coords or computational coords

    //compute delt
    double delt = time/time_slice;

    //set omega2 (rad/s)
    double omega2 = 2.0*M_PI;

    //convert degrees to radians
    alphamax = (alphamax*M_PI)/180.0;

    //determine cg (could have user input, sine 1/4 cord length for symmetric airfoil, but this is nice too)
    //however, we only want to do this at the outset, since we want our cg to be the axis of rotation!
    //create points to hold high and low extents of domains
    //set hi and low oppositely large to assure proper action of MAX and MIN
    Point lo = Point(1.0e20,1.0e20);
    Point hi = Point(-1.0e20,-1.0e20);

    int bd; //this is the generic boundary for our switch

    //determine hi and lo points of bd1 and bd2
    for (i = 0; i < 2; i++)
      {
      switch (i)
	{
	case 0: 
	  bd=bd1;
	  break;
	case 1: 
	  bd=bd2; 
	  break;
	default: 
	  printf("\nYour index is being clobbered in the switch checking boundry extents in driver.cpp!\n");
	  exit(0);
	  break;
	}      
      //check extents, first for x, then for y
      for (j = 0; j < nbs[bd]; j++)
	{
	// check both ends of each segment 
	for (k = 0; k < 2; k++)
	  {
	  //for winslow, we use computational domain boundaries
	  if (smooth_type == 0)
	    {
	    lo[0] = MIN(lo[0],nodec[bs[bd][j][k]][0]);
	    hi[0] = MAX(hi[0],nodec[bs[bd][j][k]][0]);
	    lo[1] = MIN(lo[1],nodec[bs[bd][j][k]][1]);
	    hi[1] = MAX(hi[1],nodec[bs[bd][j][k]][1]);
	    }
	  //for L-E and opt-based, we use physical domain boundaries
	  if (smooth_type == 1 || smooth_type == 2)
	    {
	    lo[0] = MIN(lo[0],nodep[bs[bd][j][k]][0]);
	    hi[0] = MAX(hi[0],nodep[bs[bd][j][k]][0]);
	    lo[1] = MIN(lo[1],nodep[bs[bd][j][k]][1]);
	    hi[1] = MAX(hi[1],nodep[bs[bd][j][k]][1]);
	    }
	  }
	}
      }

    //now, take the average of high and low y-coords and one-quarter of the length between x-coords and add to the low x-coord
    //init the center of gravity
    Point cg = Point (0.0,0.0);
    //find abs of distance between two x-vals, since they could be negative
    cg[0] = 0.25*(fabs(hi[0] - lo[0])) + lo[0];
    //check for precision here, since usually zero
    if (fabs((hi[1] + lo[1])/2.0) <= 1.0e-5)
      cg[1] = 0.0;
    else
      cg[1] = (hi[1] + lo[1])/2.0;

    //debug...print cg
    printf("\nCG = ( %lf , %lf )\n",cg[0],cg[1]);

    //again, we only want to tag bd points once, since they will remain the same node numbers throughout the oscillations
    //create special tag array to tag points on the boundaries to be moved
    int *tag2;
    tag2 = (int*)calloc(nn,sizeof(int));

    //FOR L-E: now, determine the new boundary which will be put in place of the computational nodes, which are identical to the physical nodes elsewhere and are just placeholders so we can find the deltas
    //FOR Winslow: now, determine the new boundary which will be put in place of the computational nodes AND physical nodes
    //again switch between boundaries to tag all points that are on the boundary specially, so we only change them ONCE
    for (i = 0; i < 2; i++)
      {
      switch (i)
	{
	case 0: 
	  bd=bd1;
	  break;
	case 1: 
	  bd=bd2; 
	  break;
	default: 
	  printf("\nYour index is being clobbered in the switch checking boundry extents in driver.cpp!\n");
	  exit(0);
	  break;
	 }      
      //loop through each boundary's nodes
      for (j = 0; j < nbs[bd]; j++)
	{
	// tag point at both ends of each segment, even though there is some redundancy 
	for (k = 0; k < 2; k++)
	  {
	  tag2[bs[bd][j][0]] = 1;
	  tag2[bs[bd][j][1]] = 1;
	  }
	}
      }

    //now, we want to perform this smoothing time_slices number of times, so we set up a do loop
    double current_time = 0.0; //set initial current_time to 0
    int iter_file = 0; //to concatenate with file names to separate them
    int dig; //to decide how many spaces to delete before adding new iteration tag
    char extension[bdim]; //to add iteration number to string
    do
      {
      current_time += delt; //up current time to keep up with oscillations for do loop and for delalpha formula
      iter_file++; //increment to get new filename

      double delalpha; //this will be the change in alpha from iteration to iteration
	
      //notice, when computing delalpha, by convention current_time is t+delt, so to get t we use current_time - delt
      //also, to get movement in one shot, we need to ignore periodicity
      if (time_slice == 1)
        delalpha = alphamax;
      else
        delalpha = alphamax*(sin(omega2*(current_time)) - sin(omega2*(current_time - delt)));

      //now, go through tagged nodes and change the nodec/nodep to perturbed boundary
      Vector r; //will hold Vector from cg to old nodec
      double theta; //will hold angle from midline of airfoil to old nodec with vertex at cg
      double rmag; //will hold magnitude of r
    
      for (i = 0; i < nn; i++)
	{
	if (tag2[i] == 0)
	  continue; //non boundary node
	//we move nodep for winslow smoothing OR optimization based
	if (smooth_type == 0 || smooth_type == 2)
	  {
	  r = Vector(nodep[i][0] - cg[0], nodep[i][1] - cg[1]); //find r to be used in theta and node changing
	  theta = atan2(r[1],r[0]); //find angle using atan2 to get right quadrant angle
	  rmag = r.magnitude(); //find mag of r
	  //reset nodep
	  nodep[i][0] = cg[0] + rmag*cos(theta+delalpha);
	  nodep[i][1] = cg[1] + rmag*sin(theta+delalpha);
	  }
	//we move nodec for L-E smoothing, since we are looking for a delta in our routine and the "new mesh" passed into nodec and nodep each time
	if (smooth_type == 1)
	  {
	  r = Vector(nodec[i][0] - cg[0], nodec[i][1] - cg[1]); //find r to be used in theta and node changing
	  theta = atan2(r[1],r[0]); //find angle using atan2 to get right quadrant angle
	  rmag = r.magnitude(); //find mag of r
	  //reset nodep
	  nodec[i][0] = cg[0] + rmag*cos(theta+delalpha);
	  nodec[i][1] = cg[1] + rmag*sin(theta+delalpha);
	  }
	}
      
      //these will now feed into the linear elastic, winslow, or opt based routine, which pass out new physical coords, do mesh file and create plot
      if (smooth_type == 2)
	{
	//now, pass in tag, node-cell hash, node-node hash, physical nodes and computational nodes (u,v) (second only as storage location for perturbations), tri array, and dimensions into optimization-based smoothing procedure, also we pass in cvg mostly as a token
	//also, pass iterations at which to cut off
	opt_smooth(nn, nodep, nodec, hash, NChash, nt, tri, iterin, tag, cvg);
	}

      if (smooth_type == 1)
	{
	//now, pass in mat, tag, compressed row storage arrays, physical/computational nodes (second only as storage location for perturbed boundaries), tri array, and dimensions into linear-elastic smoothing procedure, where young's modulus computer (inline procedure), weights applied (matrix build procedure), smoothing done (linear solve procedure)
	//also, pass cvg value, omega (relaxation factor), iterations at which to cut off if cvg not reached (iterin), and poisson's ratio (held constant)
	linear_elastic(nn, nodep, nodec, mdim, mat, ia, ja, iau, nt, tri, iterin, cvg, omega, tag, poisson);

	//now, we reset nodec as nodep, since we have moved the mesh and we need the proper delta to create u and v for L-E and nodec is only a placeholder
	for (i = 0; i < nn; i++)
	  {
	  nodec[i][0] = nodep[i][0];
	  nodec[i][1] = nodep[i][1];
	  }
	}

      if (smooth_type == 0)
	{
	//now, pass in mat, tag, compressed row storage arrays, physical/computational nodes, tri array, and dimensions into winslow smoothing procedure, where weights applied (matrix build procedure), smoothing done (linear solve procedure)
	//also, pass cvg value, omega (relaxation factor), iterations between recompute (iterin), iterations at which to cut off if cvg not reached (iterout)
	//poisson passed in but not used to save writing an overloaded version of mat_build
        //pass in A, B, C (forcing function coeff), whichff to tell which forcing function to create, and q to create forcing functions
        //in case we are pasing in no variables and we want no forcing functions (C = 0)
        if (nv == 0)
          {
          double **q; 
          q = (double**)calloc(nn,sizeof(double*));
          for (i = 0; i < nn; i++)
            q[i] = (double*)calloc(4,sizeof(double));
          }

	winslow(nn, nodep, nodec, mdim, mat, ia, ja, iau, nt, tri, iterin, iterout, cvg, omega, tag, poisson, A, B, C, whichff, q);

        //set initial min/max low/high for loop
        double minjac = 1.0e20;
        double maxjac = -1.0e20;
        //set counters for pos and neg jac to zero
        int posjac = 0;
        int negjac = 0;
        //holder for jacobian
        double avgjac = 0.0; 

        //check jacobians to see if smoothing successful
        for (i = 0; i < nt; i++)
          {
          //grab indices
          int t0 = tri[i][0];
          int t1 = tri[i][1];
          int t2 = tri[i][2];

          //make points
          Point pt0 = Point(nodep[t0][0],nodep[t0][1]);
          Point pt1 = Point(nodep[t1][0],nodep[t1][1]);
          Point pt2 = Point(nodep[t2][0],nodep[t2][1]);

          //make vectors
          Vector sd1 = Vector(pt1,pt2); 
          Vector sd2 = Vector(pt1,pt0);
         
          //compute jacobian
          avgjac = sd1 % sd2;

          //find min/max
          minjac = MIN(minjac,avgjac);
          maxjac = MAX(maxjac,avgjac);

          //count pos/neg jacobians
          if (avgjac >=0)
            posjac++;
          else
            negjac++;
          }

        //print results
        printf("\nMax jacobian = %lf\n",maxjac);
        printf("\nMin jacobian = %lf\n",minjac);
        printf("\nPositive jacobians = %d\n",posjac);
        printf("\nNegative jacobians = %d\n",negjac);

        //now, free placeholder mem
        if (nv == 0)
          {
          for (i = 0; i < nn; i++)
            free(q[i]);
          free(q);
          }
	}

      //now, write out new physical mesh file 
      //take off _%d%d.dat from old file name and reset to '\0' to re-concatenate _%d%d+1.mesh
      //find starting position of extension
      i = strlen(filename) - strlen(strstr(filename,".dat"));
      //set dig to 3 spaces since we will add zeros to single digit numbers
      dig = 3;
      //reset to '\0' dig further back to get rid of number and _
      filename[i-dig] = '\0';
      //reconcatenate
      extension[0] = '\0'; //reset extension
      dig = 0; //reset dig to zero for another use
      if (iter_file/10 == 0) //doing integer math...if less than 10 iterations, we add 0 in front of iter_file
	dig = 1;
      //now, if dig, we add 0 in front of iter_file
      if (dig)
	sprintf(extension,"_0%d",iter_file);  //we put iteration in extension
      else
	sprintf(extension,"_%d",iter_file);  //we put iteration in extension
      strcat(filename,extension); //add extension
      strcat(filename,".mesh"); //add file extension
      printf("\nOutput Mesh File filename = <%s>\n",filename);
      
      //now, open file for write
      if ((fp = fopen(filename,"w")) == 0)
	{
	printf("\nError opening file <%s>.",filename);
	exit(0);
	}

      //now, fill in each section with current physical info (noting that connectivity from computational mesh still intact)
      fprintf(fp,"#Number of grid points\n");
      fprintf(fp,"%d\n",nn);
      for (i = 0; i < nn; i++)
	fprintf(fp,"%16.10e %16.10e\n",nodep[i][0],nodep[i][1]);

      fprintf(fp,"#Number of blocks\n");
      fprintf(fp,"%d\n",nblk);
      
      //note, tri will need to be re-fortran indexed
      fprintf(fp,"#Number of triangles\n");
      fprintf(fp,"%d\n",nt);
      for (i = 0; i < nt; i++)
	fprintf(fp,"%d %d %d\n",tri[i][0]+1,tri[i][1]+1,tri[i][2]+1);
    
      //note, quad will need to be re-fortran indexed
      fprintf(fp,"#Number of quads\n");
      fprintf(fp,"%d\n",nq);
      if (nq > 0)
	{
	for (i = 0; i < nq; i++)
	  fprintf(fp,"%d %d %d %d\n",quad[i][0]+1,quad[i][1]+1,quad[i][2]+1,quad[i][3]+1);
	}

      fprintf(fp,"#Number of boundaries\n");
      fprintf(fp,"%d\n",nb);

      //note, bs will need to be re-fortran indexed
      for (i = 0; i < nb; i++)
	{
	fprintf(fp,"#Number of edges on boundary %d\n",i);
	fprintf(fp,"%d\n",nbs[i]);
	for (j = 0; j < nbs[i]; j++)
	  {
	  fprintf(fp,"%d %d\n",bs[i][j][0]+1,bs[i][j][1]+1);
	  }
	}

      //I WANT SOLVER TO BEGIN FROM SCRATCH!

      // write number of constants
      fprintf(fp,"#Number of constants\n");
      fprintf(fp,"%d\n",nc);

      if (nc > 0)
        {  
        // write out constants (we will need to know what each is in order, a priori)
        for (i = 0; i < nc; i++)
          fprintf(fp,"%lf\n",constants[i]);
        }

      // write out number of vars (usually four, but we could code for more if needed)
      fprintf(fp,"#Number of variables\n");
      fprintf(fp,"%d\n",0);

      /*//we do not have it set up to read in variable names, so we just write out what we know they are a priori
      fprintf(fp,"density\n");
      fprintf(fp,"x-momentum\n");
      fprintf(fp,"y-momentum\n");
      fprintf(fp,"total energy\n");

      //now, write out q....hard wired to four vars  
      if (nv > 0)
        {  
        // write out q
        for (i = 0; i < nn; i++)
          fprintf(fp,"%16.10e %16.10e %16.10e %16.10e\n",q[i][0],q[i][1],q[i][2],q[i][3]);
        }*/

      fclose(fp);

      //finally, write GNUPLOT file
      //take off .mesh from old file name and reset to '\0' to re-concatenate .dat
      //find starting position of file extension
      i = strlen(filename) - strlen(strstr(filename,".mesh"));
      //reset to '\0'
      filename[i] = '\0';
      //reconcatenate
      strcat(filename,".dat");
      printf("\nOutput Plot File Filename = <%s>\n",filename);
      // Open file for write
      if ((fp = fopen(filename,"w")) == 0)
	{
	printf("\nError opening file <%s>.",filename);
	exit(0);
	}
      for (i=0; i < nt; i++)
	{
	n0 = tri[i][0];
	n1 = tri[i][1];
	n2 = tri[i][2];
	fprintf(fp,"%19.10e %19.10e 0.0\n",  nodep[n0][0],nodep[n0][1]);
	fprintf(fp,"%19.10e %19.10e 0.0\n",  nodep[n1][0],nodep[n1][1]);
	fprintf(fp,"%19.10e %19.10e 0.0\n",  nodep[n2][0],nodep[n2][1]);
	fprintf(fp,"%19.10e %19.10e 0.0\n\n",nodep[n0][0],nodep[n0][1]);
	}
    
      fclose(fp);
   
      } while (current_time < time); //stop once we have done prescribed number of time steps


    //finally, let's compute RMSerror between each node of orig and final meshes
    double RMSerr = 0.0;
    for (i = 0; i < nn; i++)
      RMSerr += (nodep[i][0] - tnode[i][0])*(nodep[i][0] - tnode[i][0]) + (nodep[i][1] - tnode[i][1])*(nodep[i][1] - tnode[i][1]);
    //compute final RMS remembering that we have 2*nn elements since we have errors in x and y
    RMSerr = sqrt(RMSerr/(2*nn));

    //print to screen
    printf("\nRMSerr between original and final meshes is %16.10e\n",RMSerr);
    
    //free memory, reset to null pointers (for this if block)
    if (smooth_type == 0 || smooth_type == 1)
      {
      freenull(ia);
      freenull(ja);
      freenull(iau);
      freenull(mat);
      }
    freenull(tag);
    freenull(tag2);
    freenull(tnode);
    for (i = 0; i < nn; i++)
      delete NChash[i];
    delete[] NChash;
    if (smooth_type == 2)
      {
      //this is done earlier in l-e and winslow, but not until end in opt-based
      for (i = 0; i < nn; i++)
        delete hash[i];
      delete[] hash;
      }
    }

  if (operation == 1)
    {
    //read in number of functions
    int test_funct = 0;
    printf("\nEnter number of functions to be tested for refinement ->");
    fgets(buff,bdim,stdin); //read in test_funct
    sscanf(buff,"%d",&test_funct); //place in variable
    //error checking
    if (test_funct <= 0 || (test_funct > 1 && nv ==0) || test_funct == 2 || test_funct > 3)
      {
      printf("\nYou cannot refine based on %d functions with %d variables.  Exiting....", test_funct, nv);
      fflush(stdout);
      exit(0);
      }

    int ft = 0; //which function to use var toggle
    //read in which analytic function is to be used
    if (nv == 0)
      {
      printf("\nEnter 0 to use the circle analytic function or 1 to use the cosine analytic function ->");
      fgets(buff,bdim,stdin); //read in ft
      sscanf(buff,"%d",&ft); //place in variable
      //error checking
      if (ft < 0 || ft > 1)
        {
        printf("\nYou cannot refine based on analytic function number %d.  Exiting....", ft);
        fflush(stdout);
        exit(0);
        }
      }
    //read in which solution based function(s) is/are to be used
    if (nv == 4)
      {
      printf("\nEnter 2 to use the pressure function, 3 to use the velocity magnitude function, 4 to use the mach number function, or 5 to use all three ->");
      fgets(buff,bdim,stdin); //read in ft
      sscanf(buff,"%d",&ft); //place in variable
      //error checking
      if (ft < 2 || ft > 5 || (ft == 5 && test_funct != 3))
        {
        printf("\nYou cannot refine based on function number %d with %d total functions.  Exiting....", ft,test_funct);
        fflush(stdout);
        exit(0);
        }
      }


    //read in type of refinement funct
    int funct_type = 0;
    printf("\nEnter 0 for gradient-based or 1 for edge-based refinement ->");
    fgets(buff,bdim,stdin); //read in funct_type
    sscanf(buff,"%d",&funct_type); //place in variable
    //error checking
    if (funct_type < 0 || funct_type > 1)
      {
      printf("\nYou must choose 0 or 1 for gradient-based or edge-based refinement.  Exiting....");
      fflush(stdout);
      exit(0);
      }

    //read in num of std deviations needed to refine
    double Cr[3];
    for (i = 0; i < test_funct; i++)
      {
      printf("\nEnter number of standard deviations needed for refining function %d ->",i);
      fgets(buff,bdim,stdin); //read in Cr
      sscanf(buff,"%lf",&Cr[i]); //place in variable
      //error checking
      if (Cr[i] <= 0.0)
        {
        printf("\nYou cannot use %lf std deviations for refining.  Exiting....",Cr[i]);
        fflush(stdout);
        exit(0);
        }
      }

    //read in num of std deviations needed to coarsen
    double Cc[3];
    for (i = 0; i < test_funct; i++)
      {
      printf("\nEnter number of standard deviations needed for coarsening function %d ->",i);
      fgets(buff,bdim,stdin); //read in Cc
      sscanf(buff,"%lf",&Cc[i]); //place in variable
      //error checking
      if (Cc[i] <= 0.0)
        {
        printf("\nYou cannot use %lf std deviations for coarsening.  Exiting....",Cc[i]);
        fflush(stdout);
        exit(0);
        }
      } 

    //read in power to raise length to for refinement function
    double p[3];
    for (i = 0; i < test_funct; i++)
      {
      printf("\nEnter power to raise length to for refinement function %d ->",i);
      fgets(buff,bdim,stdin); //read in p
      sscanf(buff,"%lf",&p[i]); //place in variable
      //error checking
      if (p[i] < 1.0)
        {
        printf("\nYou cannot use %lf as the power in computing the refinement function.  Exiting....",p[i]);
        fflush(stdout);
        exit(0);
        }
      } 
      
    //read in whether or not to retriangulate
    remesh = 0;
    printf("\nEnter 1 to retriangulate or 0 to leave as is (assuming no coarsening done) ->");
    fgets(buff,bdim,stdin); //read in remesh
    sscanf(buff,"%d",&remesh); //place in variable
    //error checking
    if (remesh < 0 || remesh > 1)
      {
      printf("\nYou must choose 0 for no retriangulation or 1 for a new Dealuney mesh.  Exiting....",p);
      fflush(stdout);
      exit(0);
      }

    //we must expand tri to higher order connectivity
    for (i = 0; i < nt; i++)
      {
      tri[i] = (int*)realloc((void*)tri[i],6*sizeof(int));
      //init all entries to -1
      tri[i][3] = tri[i][4] = tri[i][5] = -1;
      }

    //allocate for node to node hash to be passed into routine and filled in
    List **hash;
    hash = new List*[nn]; //each node has a list of nodes it's attached to
    for (i=0; i < nn; i++)
      hash[i] = new List(); //make each node's list

    //now, create hash table to be passed to compressed row storage routine
    Lhash(nn, nt, tri, hash);

    //now, we build a nabor connectivity
    //init Linked_Lists
    Linked_List **Nhash;
    Linked_List *mknbr;

    //first, initialize node to element hash table
    Nhash = new Linked_List*[nn];

    for (i = 0; i < nn; i++)
      Nhash[i] = new Linked_List();

    //next, initialize list of all triangles (to keep track of nabor making needs in make_nbr)
    mknbr = new Linked_List();

    //now, insert each tri into mknbr and simultaneously create node to element hash
    for (i = 0; i < nt; i++)
      {
      mknbr->Insert(i); //insert tri into list
      for (j = 0; j < 3; j++)
        {
        n = tri[i][j]; //use n as placeholder for current node index (loop over all three)
        Nhash[n]->Insert(i);
        }
      }

    //now, create nabor connectivity to be passed in init to -1 and out full
    int **nbr;
    nbr = (int**)calloc(nt,sizeof(int*));
    for (i = 0; i < nt; i++)
      nbr[i] = (int*)calloc(3,sizeof(int));
    for (i = 0; i < nt; i++)
      nbr[i][0]=nbr[i][1]=nbr[i][2]=-1; //init to -1

    //now, we have all our hash tables and we can create nabor connectivity
    make_nbrs_HO(mknbr, nn, nt, tri, nbr, Nhash);

    //now, we can delete mknbr, since only used for convenience of keeping track in make_nbr
    delete mknbr;

    //create tag array of boundry info, no need to send nb, nbs, bs to subroutines
    int *tag = (int*)calloc(nn,sizeof(int));
 
    //init tag to ones, we need tagged boundaries for both gradients and areas since we deal with closing the contours there differently
    for (i = 0; i < nn; i++)
      tag[i] = 1;

    //1 for interior, 0 for boundary
    for (i = 0; i < nb; i++)
      {
      for (j = 0; j < nbs[i]; j++)
        {
        for (k = 0; k < 2; k++)
          {
          if (tag[bs[i][j][k]] == 1)
            tag[bs[i][j][k]] = 0; //just to be sure we don't waste time overwriting
          }
        }
      }

    //also, before doing coarsening, we define map array to help with deletion
    //we also do this here due to the need to realloc
    int *map = (int*)calloc(nn,sizeof(int));

    //init to -1
    for (i = 0; i < nn; i++)
      map[i] = -1;

    //now, take all boundary nodes and set to 0 (keep)
    for (i = 0; i < nn; i++)
      if (tag[i] == 0)
        map[i] = 0;
  
    //now, we pass info into subdivision routine to compute functions, gradients, determine refinement status
    int newnodes = subdivision(q, nv, test_funct, tri, nt, nodep, nn, funct_type, Cr, Cc, p, hash, nbr, tag, map, ft);

    //we will eventually use this list to remap to new global node numbers and send only these to trimesh to re-triangulate
    //now, we set map of new nodes to zero and realloc map
    //also, if there are nodes deleted, we will do that later
    map = (int*)realloc((void*)map,newnodes*sizeof(int));
    for (i = nn; i < newnodes; i++)
      map[i] = 0; //we always keep refined nodes
  
    //now, we want to generate new nodes and realloc nodep
    //again, we will get rid of deleted nodes before we send to trimesh
    nodep = (double**)realloc((void*)nodep,newnodes*sizeof(double*));
    for (i = nn; i < newnodes; i++)
      nodep[i] = (double*)calloc(2,sizeof(double));

    //now, we loop through triangles and store physical coords of new nodes in nodep
    //first, declare some placeholder vars
    int nd0, nd1; //node index placeholders
    Point p0, p1, p2, pt0, pt1, newpt; //point class placeholders
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

      //now, switch over all three sides to check for new nodes
      for (j = 0; j < 3; j++)
        {
        switch(j)
          {
          case 0:
            pt0 = p0;
            pt1 = p1;  
          break;
          case 1:
            pt0 = p1;
            pt1 = p2; 
          break;
          case 2:
            pt0 = p2;
            pt1 = p0;
          break;
          default:
            printf("\nYou have a triangle with more than three sides according to the loop counter in driver.cpp!\n");
            fflush(stdout);
            exit(0);
          break;
          }
        //now, define new point
        if ((n = tri[i][j+3]) >= 0)
          {
          newpt = (pt0 + pt1)/2.0;
          nodep[n][0] = newpt[0];
          nodep[n][1] = newpt[1];
          }
        }
      }

    //now, transfer solution by extending q and using averaging to find q at new points
    if (nv > 0)
      {
      q = (double**)realloc((void*)q,newnodes*sizeof(double*));
      for (i = nn; i < newnodes; i++)
        q[i] = (double*)calloc(nv,sizeof(double));
      
      //we will use triangle conn to find new q val
      for (i = 0; i < nt; i++)
        {
        //use node indices for cleanliness
        n0 = tri[i][0];
        n1 = tri[i][1];
        n2 = tri[i][2];
      
        //now, switch over all three sides to update q at any new nodes
        for (j = 0; j < 3; j++)
          {
          switch(j)
            {
            case 0:
              nd0 = n0;
              nd1 = n1;  
            break;
            case 1:
              nd0 = n1;
              nd1 = n2; 
            break;
            case 2:
              nd0 = n2;
              nd1 = n0;
            break;
            default:
              printf("\nYou have a triangle with more than three sides according to the loop counter in driver.cpp!\n");
              fflush(stdout);
              exit(0);
            break;
            }
          //now, define new q vals
          if ((n = tri[i][j+3]) >= 0)
            for (k = 0; k < nv; k++)
              q[n][k] = (q[nd0][k] + q[nd1][k])/2.0;
          }
        }
      }

    //now, reset nn to newnodes and old_nn to orig nn for memory freeing
    old_nn = nn;
    nn = newnodes;

    //now, we need to update boundaries by looping thru bd seg, using node_cell hash to find tri, looking for new node on edge, and reallocing as we go
    //linked nodes to read thru linked list (node_cell hash)
    Linked_Node *hd0, *hd1;
    for (i = 0; i < nb; i++)
      {
      for (j = 0; j < nbs[i]; j++)
        {
        //set nodes to look for
        n0 = bs[i][j][0];
        n1 = bs[i][j][1];

        //we don't need to look thru bs with newly created nodes, since they cannot be split
        if (n0 >= old_nn || n1 >= old_nn)
          continue;
         
        //init m (triangle with bs in it) to -1
        m = -1;
        
        // each element in n0's list is a triangle with that node in it
        hd0 = Nhash[n0]->head;

        //while there are elements that have n0 in them and m not set, keep looking
        while (hd0 && m < 0)
          {
          //cell from n0 Nhash
          n = hd0->data;

          //now, look to other node on bd seg's Nhash to see if it is also in element n
          hd1 = Nhash[n1]->head;

          //while there are elements in n1's Nhash and m unset
          while (hd1 && m < 0)
            {
            //if the element in n1's hash matches the element in n0's hash, we have the element that constains bd seg
            //set triangle that has both bs nodes as m
            if (hd1->data == n)
              m = n;

            //move to the next element in n1's hash, even if unneeded
            hd1 = hd1->next;
            }
          //move to the next element in n0's hash, even if unneeded
          hd0 = hd0->next;
          }

        //check that a triangle has been found...if not, we have a problem!
        if (m < 0)
          {
          printf("\nBoundary segment %d on boundary %d is NOT in a triangle!  Exiting....\n",j,i);
          fflush(stdout);
          exit(0);
          }

        //now that we have a triangle that bs is in, look to see if we have a new node on the given side
        //run thru a switch
        for (s = 0; s < 3; s++)
          {
          switch(s)
            {
            case 0:
              nd0 = tri[m][0];
              nd1 = tri[m][1];  
            break;
            case 1:
              nd0 = tri[m][1];
              nd1 = tri[m][2]; 
            break;
            case 2:
              nd0 = tri[m][2];
              nd1 = tri[m][0];
            break;
            default:
              printf("\nYou have a triangle with more than three sides according to the loop counter in driver.cpp!\n");
              fflush(stdout);
              exit(0);
            break;
            }
          //look thru sides, if find side (with new point) then reset bd seg, realloc bs
          if (nd0 == n0 && nd1 == n1)
            {
            if (tri[m][s+3] >= 0)
              {
              nbs[i]++; //increment number of boundaries
              bs[i] = (int**)realloc((void*)bs[i],nbs[i]*sizeof(int*)); //realloc number of bs for bd i          
              for (k = nbs[i]-1; k < nbs[i]; k++)
                bs[i][k] = (int*)calloc(2,sizeof(int)); //add spots for two more indices at end of list
              //now, place first half in orig spot and second half in new spot
              //be careful to keep ordering straight
              bs[i][j][0] = n0;
              bs[i][j][1] = tri[m][s+3];
              bs[i][nbs[i]-1][0] = tri[m][s+3];
              bs[i][nbs[i]-1][1] = n1;
              }
            }
          if (nd0 == n1 && nd1 == n0)
            {
            printf("\nBoundary segment %d on boundary %d is out of order with triangle side %d on triangle &d!\n",j,i,s,m);
            fflush(stdout);
            exit(0);
            }
          }
        }
      }
    
    //we can free nbr now since we no longer need it, and nt and old_nt will change
    for (i = 0; i < nt; i++)
      free(nbr[i]);
    freenull(nbr);

    //now, we count the number of new triangles created so as to realloc tri
    //first, we need to set old_nt to current nt (for loop), since nt will be incremented as we go
    int old_nt = nt;

    //printf("\nold_nt = %d",old_nt);
    //fflush(stdout);

    //also, let's create a triangle array (original style connectivity) which is hardwired to 4 triangles (most we can get in 2-D)
    const int cdim = 4;
    int tri_piece[cdim][3]; //no need to free since static alloc
    //also, let's create a temp array to hold our current triangles conn
    int conn[6]; //no need to free since static alloc
    
    //even if we are going to retriangulate, we do this so that we have a better idea of how many triangles we may have

    for (i = 0; i < old_nt; i++)
      {
      //first, ready the current tri higher order conn to be passed in
      for (j = 0; j < 6; j++)
        conn[j] = tri[i][j];

      //now, reset the tri_piece array
      for (j = 0; j < cdim; j++)
        for (k = 0; k < 3; k++)
          tri_piece[j][k] = 0;
      
      //now, pass into refine tri only to determine new nt (use l as dummy, since if returns 1, no new triangles!)
      l = refine_tri(conn, cdim, tri_piece, nodep);

      //now, if l > 1, add l - 1 to tri, since we will be deleting the mother triangle each time
      if (l > 1)
        nt += (l-1);
      }

    //printf("\nnt = %d",nt);
    //fflush(stdout);

    //now, we realloc tri
    tri = (int**)realloc((void*)tri,nt*sizeof(int*));
    for (i = old_nt; i < nt; i++)
      tri[i] = (int*)calloc(6,sizeof(int));

    //also reset new tri entries 3, 4, 5 to -1
    for (i = old_nt; i < nt; i++)
      for (j = 3; j < 6; j++)
        tri[i][j] = -1;

    //we set a flagC to find out if coarsening must occur
    flagC = 0;
    //now check for coarsening
    for (i = 0; i < nn; i++)
      {
      if (map[i] == -1)
        flagC++; //will reuse for realloc      
      }
    //we only do this if remesh set to 0 AND we have NOT coarsened
    if (remesh == 0 && !flagC)
      {
      //now, we add in new triangles
      for (i = 0; i < old_nt; i++)
        {
        //first, ready the current tri higher order conn to be passed in
        for (j = 0; j < 6; j++)
          conn[j] = tri[i][j];

        //now, reset the tri_piece array
        for (j = 0; j < cdim; j++)
          for (k = 0; k < 3; k++)
            tri_piece[j][k] = 0;
      
        //now, pass into refine tri to determine new connectivity (refine_tri gives global nodes)
        //l is a dummy var that holds new triangles gen to help loop counter and conditional for adding new tri
        l = refine_tri(conn, cdim, tri_piece, nodep);

        //now, take tri_piece array and add first triangle into old tri spot, and the rest to the bottom of the list
        if (l > 1)
          {
          //first new triangle goes in old tri spot
          for (j = 0; j < 3; j++)
            tri[i][j] = tri_piece[0][j];
          
          //now, add the rest of the new triangles to tri starting at old_nt
          for (j = 1; j < l; j++)
            {
            for (k = 0; k < 3; k++)
              {
              //reset each element in new tri conn to tri_piece passed out from refine_tri
              tri[old_nt][k] = tri_piece[j][k];
              }
            old_nt++; //increment old_nt AFTER adding EACH triangle
            }
          }
        } 
      }

    //finally, reduce conn to lower order no matter what!
    for (i = 0; i < nt; i++)
      tri[i] = (int*)realloc((void*)tri[i],3*sizeof(int));
     
    // we only bother to remap if we had deletions   
    if (flagC)
      {      
      //now, we wish to remap the undeleted and refined and boundary nodes
      int nodemap = 0; //start numbering at zero
      for (i = 0; i < nn; i++)
        {
        if (map[i] == 0)
          {
          map[i] = nodemap;
          nodemap++;
	  }
        }
      
      //now, we need to remap the physical nodes and realloc
      for (i = 0; i < nn; i++)
        {
        if (map[i] >= 0 && map[i] != i)
          {   
          nodep[map[i]][0] = nodep[i][0]; //move x node back
	  nodep[map[i]][1] = nodep[i][1]; //move y node back
	  }
        }
	
      //printf("\nnn = %d\n",nn);	
      //printf("\nflagC = %d\n",flagC);
	
      //now realloc nodep
      for (i = (nn-flagC); i < nn; i++)
        freenull(nodep[i]);
      nodep = (double**)realloc((void*)nodep,(nn-flagC)*sizeof(double*));
    
      //now, we need to remap bs
      for (i = 0; i < nb; i++)
        for (j = 0; j < nbs[i]; j++)
          for (k = 0; k < 2; k++)
	    bs[i][j][k] = map[bs[i][j][k]];
	    
      //now, we need to remap q and realloc
      for (i = 0; i < nn; i++)
        {
	for (j = 0; j < nv; j++)
	  {
          if (map[i] >= 0 && map[i] != i)
            {   
            q[map[i]][j] = q[i][j]; //move q_j node back
	    }
	  }
        }
	
      //now realloc q
      for (i = (nn-flagC); i < nn; i++)
        freenull(q[i]);
      q = (double**)realloc((void*)q,(nn-flagC)*sizeof(double*));
      
      //reset nn
      nn -= flagC;
      //printf("\nnew nn = %d\n",nn);
      }
	  
    //now, we wish to re trimesh if we had deletions or remesh == 1
    //note that the trimesh that is passed in is void and will be redone completely and reallocated if necessary
    //however, if less triangles are generated, than all entries past nt are garbage, but it is the end of the program and we needn't realloc!
    if (remesh == 1 || flagC)
      {
      //now, we need to place nodep into temp arrays for passing to trimesh
      double *x, *y;
      x = (double*)calloc(nn,sizeof(double));
      y = (double*)calloc(nn,sizeof(double));
      
      for (i = 0; i < nn; i++)
        {
	x[i] = nodep[i][0];
	y[i] = nodep[i][1];
	}
      
      // iterative loop to allocate space for triangle indices if more is needed
      tdim = -1; //start at -1 to assure goes thru loop
      // if we run out of space, trimesh returns -1 and asks to realloc by adding increments on nn
      do
        {
        tdim = trimesh(nn, nt, nb, nbs, bs, x, y, tri);
 
        if (tdim < 0)
          {
          printf("\nExpanding triangle dimension from %d to %d",nt,nt+nn);
          fflush(stdout);
        
          // increment tri array dimension by number of nodes
          nt += nn;
        
          //now, we realloc tri..even though the routine will start from scratch
          tri = (int**)realloc((void*)tri,nt*sizeof(int*));
          for (i = nt-nn; i < nt; i++)
            tri[i] = (int*)calloc(3,sizeof(int));
          }
        } while (tdim < 0);
	
      //nt is already reset for allocation purposes, but we wish to reset it so if we have extra garbage at the end, we don't print it out
      nt = tdim;            
      }
              
    //now, write out new physical mesh file 
    strcat(filename,".mesh"); //add file extension
    printf("\nOutput Mesh File filename = <%s>\n",filename);
    
    //now, open file for write
    if ((fp = fopen(filename,"w")) == 0)
      {
      printf("\nError opening file <%s>.",filename);
      exit(0);
      }

    //now, fill in each section with current physical info (noting that connectivity from computational mesh still intact)
    fprintf(fp,"#Number of grid points\n");
    fprintf(fp,"%d\n",nn);
    for (i = 0; i < nn; i++)
      fprintf(fp,"%16.10e %16.10e\n",nodep[i][0],nodep[i][1]);

    fprintf(fp,"#Number of blocks\n");
    fprintf(fp,"%d\n",nblk);
    
    //note, tri will need to be re-fortran indexed
    fprintf(fp,"#Number of triangles\n");
    fprintf(fp,"%d\n",nt);
    for (i = 0; i < nt; i++)
      fprintf(fp,"%d %d %d\n",tri[i][0]+1,tri[i][1]+1,tri[i][2]+1);
  
    //note, quad will need to be re-fortran indexed
    fprintf(fp,"#Number of quads\n");
    fprintf(fp,"%d\n",nq);
    if (nq > 0)
      {
      for (i = 0; i < nq; i++)
	fprintf(fp,"%d %d %d %d\n",quad[i][0]+1,quad[i][1]+1,quad[i][2]+1,quad[i][3]+1);
      }

    fprintf(fp,"#Number of boundaries\n");
    fprintf(fp,"%d\n",nb);

    //note, bs will need to be re-fortran indexed
    for (i = 0; i < nb; i++)
      {
      fprintf(fp,"#Number of edges on boundary %d\n",i);
      fprintf(fp,"%d\n",nbs[i]);
      for (j = 0; j < nbs[i]; j++)
	{
	fprintf(fp,"%d %d\n",bs[i][j][0]+1,bs[i][j][1]+1);
	}
      }

    // write number of constants
    fprintf(fp,"#Number of constants\n");
    fprintf(fp,"%d\n",nc);

    if (nc > 0)
      {  
      // write out constants (we will need to know what each is in order, a priori)
      for (i = 0; i < nc; i++)
        fprintf(fp,"%lf\n",constants[i]);
      }

    // write out number of vars (usually four, but we could code for more if needed)
    fprintf(fp,"#Number of variables\n");
    fprintf(fp,"%d\n",nv);

    //we do not have it set up to read in variable names, so we just write out what we know they are a priori
    fprintf(fp,"density\n");
    fprintf(fp,"x-momentum\n");
    fprintf(fp,"y-momentum\n");
    fprintf(fp,"total energy\n");

    //now, write out q....hard wired to four vars  
    if (nv > 0)
      {  
      // write out q
      for (i = 0; i < nn; i++)
        fprintf(fp,"%16.10e %16.10e %16.10e %16.10e\n",q[i][0],q[i][1],q[i][2],q[i][3]);
      }

    fclose(fp);

    //finally, write GNUPLOT file
    //take off .mesh from old file name and reset to '\0' to re-concatenate .dat
    //find starting position of file extension
    i = strlen(filename) - strlen(strstr(filename,".mesh"));
    //reset to '\0'
    filename[i] = '\0';
    //reconcatenate
    strcat(filename,".dat");
    printf("\nOutput Plot File Filename = <%s>\n",filename);
    // Open file for write
    if ((fp = fopen(filename,"w")) == 0)
      {
      printf("\nError opening file <%s>.",filename);
      exit(0);
      }
    for (i = 0; i < nt; i++)
      {
      n0 = tri[i][0];
      n1 = tri[i][1];
      n2 = tri[i][2];
      fprintf(fp,"%19.10e %19.10e 0.0\n",  nodep[n0][0],nodep[n0][1]);
      fprintf(fp,"%19.10e %19.10e 0.0\n",  nodep[n1][0],nodep[n1][1]);
      fprintf(fp,"%19.10e %19.10e 0.0\n",  nodep[n2][0],nodep[n2][1]);
      fprintf(fp,"%19.10e %19.10e 0.0\n\n",nodep[n0][0],nodep[n0][1]);
      }
  
    fclose(fp);

    //free up memory (for this if block)
    for (i = 0; i < old_nn; i++)
      delete Nhash[i];
    delete [] Nhash;
    for (i = 0; i < old_nn; i++)
      delete hash[i];
    delete[] hash;
    if (nv > 0)
      {
      for (i = 0; i < nn; i++)
        free(q[i]);
      freenull(q);
      }
    if (nc > 0)
      freenull(constants);
    }

  //free memory, reset to null pointers (for whole program)
  //depending on if we smoothed or refined, we need to delete a different amount of nodec
  if (operation == 0)
    {
    for (i = 0; i < nn; i++)
      free(nodec[i]);
    freenull(nodec);
    }
  if (operation == 1)
    {
    for (i = 0; i < old_nn; i++)
      free(nodec[i]);
    freenull(nodec);
    }
  for (i = 0; i < nn; i++)
    free(nodep[i]);
  freenull(nodep);
  //if we did retriangulation, we need to be sure to free the right amount of memory, so we use tdim
  if (remesh == 1 || flagC)
    {
    for (i = 0; i < tdim; i++)
      free(tri[i]);
    freenull(tri);
    }
  else
    {
    for (i = 0; i < nt; i++)
      free(tri[i]);
    freenull(tri);
    }
  if (nq > 0)
    {
    for (i = 0; i < nq; i++)
      free(quad[i]);
    freenull(quad);
    }
  for (i = 0; i < nb; i++)
    for (j = 0; j < nbs[i]; j++)
      free(bs[i][j]);
  for (i = 0; i < nb; i++)
    free(bs[i]);
  freenull(bs);
  freenull(nbs);

  return(0);
}
