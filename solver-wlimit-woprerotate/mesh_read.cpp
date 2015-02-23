#include "Mesh.h"

void Mesh::mesh_read(char filename[])
  { 
  int i, j, k; //counters
  const int bdim = 132; //buffer size
  char buff[bdim]; //buffer for reading in
  FILE *fp = 0; //file pointer
  
  //check to be sure physical mesh file can be opened
  if ((fp=fopen(filename,"r")) == 0)
    {
    printf("\nCouldn't open physical mesh file <%s>\n",filename);
    fflush(stdout);
    exit(0);
    }

  // read number of nodes
  fgets(buff,bdim,fp); //skip commment
  fgets(buff,bdim,fp); //read nn
  sscanf(buff,"%d",&nn);
  printf("\nNumber of nodes = %d",nn);
  fflush(stdout); //be sure prints to screen

  // allocate for nodes, using calloc, due to possible realloc
  nodep = (double**)calloc(nn,sizeof(double*));
  for (i = 0; i < nn; i++) 
    nodep[i] = (double*)calloc(2,sizeof(double));

  // read in coordinates
  for (i = 0; i < nn; i++)
    {
    fgets(buff,bdim,fp); //read node coords to be placed in row of node number (c-indexed)
    sscanf(buff,"%lf %lf",&(nodep[i][0]),&(nodep[i][1]));
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
      q[i] = (double*)calloc(4,sizeof(double));

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

  //we do not need further data from file, so we close it
  fclose(fp);

  //write GNUPLOT file for debugging (to be sure read in correctly)
  //get ready to reconcatenate with .dat
  i = strlen(filename) - strlen(strstr(filename,".mesh")); 
  //reset filename to '\0' before extension
  filename[i] = '\0';
  //add file extension
  strcat(filename,".dat"); 
  printf("\nOutput Plot File Name = <%s>\n",filename);

  // Open file for write
  if ((fp = fopen(filename,"w")) == 0)
    {
    printf("\nError opening file <%s>.",filename);
    fflush(stdout);
    exit(0);
    }
 
  for (i=0; i < nt; i++)
    {
    int n0 = tri[i][0];
    int n1 = tri[i][1];
    int n2 = tri[i][2];
    fprintf(fp,"%19.10e %19.10e 0.0\n",  nodep[n0][0],nodep[n0][1]);
    fprintf(fp,"%19.10e %19.10e 0.0\n",  nodep[n1][0],nodep[n1][1]);
    fprintf(fp,"%19.10e %19.10e 0.0\n",  nodep[n2][0],nodep[n2][1]);
    fprintf(fp,"%19.10e %19.10e 0.0\n\n",nodep[n0][0],nodep[n0][1]);
    }

  fclose(fp);

  //reset fp to NULL
  fp = 0;

  return;
  }

