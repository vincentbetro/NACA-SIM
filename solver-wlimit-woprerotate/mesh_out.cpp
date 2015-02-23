#include "Mesh.h"

void Mesh::mesh_out(char filename[])
  { 
  int i, j, k; //counters
  FILE *fp = 0; //file pointer

  //rename output mesh file 
  //take off .mesh from old file name and reset to '\0' to re-concatenate _out.mesh
  //find starting position of extension
  i = strlen(filename) - strlen(strstr(filename,".dat"));
  //reset to '\0'
  filename[i] = '\0';
  
  //now, write out new physical mesh file 
  strcat(filename,"_out.mesh"); //add file extension
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

  //hard wire nv=4 at this point, since we will read in with no q vals initially, but write out with them
  nv = 4;

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

  //we do not need further data from file, so we close it
  fclose(fp);

  //write GNUPLOT file for debugging (to be sure read in correctly)
  //get ready to reconcatenate with .dat
  /*i = strlen(filename) - strlen(strstr(filename,".mesh")); 
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
  fp = 0;*/

  return;
  }
