#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "Mesh.h"

//for debugging
void wait ( int seconds )
  {
  clock_t endwait;
  endwait = clock () + seconds * CLOCKS_PER_SEC ;
  while (clock() < endwait) {}
  }

int main(int argcs, char* pArgs[])
  { 
  //first, we need to initialize the flow by reading in Mach, alpha, gamma
  //declare variables to be read into
  double alpha, alpha1, xmach, gma; //angle of attack, mach number, gamma (to eventually be placed in mesh_obj)
  //declare vars to convert deg to rad
  double pi = 4.0*atan(1.0); //could use M_PI with math library, too
  double conv = 180.0/pi; //converts from degrees to radians
  //declare var to hold initial primitive variable values
  double rho0, c0, p0, u0, v0, ei0, et0; //non-dimensionalized prim var init vals
  const int bdim = 132; //buffer size
  char buff[bdim]; //buffer for reading in
  int i, j, k, l, n, n0, n1, n2; //counters, dummy indices
  double cp, p, x, y; //dummy var for printing tecplot
  //declare timing vars
  struct timeval start_time;
  struct timeval end_time;
  double delta_time;

  //read in alpha, xmach, gamma
  //myoot if using a restart file
  printf("\nEnter angle of attack (alpha) ->");
  fgets(buff,bdim,stdin); //read in angle of attack
  sscanf(buff,"%lf",&alpha1); //place in variable
  printf("\nEnter Mach ->");
  fgets(buff,bdim,stdin); //read in mach
  sscanf(buff,"%lf",&xmach); //place in variable
  printf("\nEnter gamma ->");
  fgets(buff,bdim,stdin); //read in gamma
  sscanf(buff,"%lf",&gma); //place in variables

  //since meshes will come in rotated, no need to rotate free stream
  //VCB:  note..changed for flux limiter calcs
  alpha = alpha1;

  //recast alpha1
  int alpha2 = (int) alpha1;

  //now, set initial values of primitive vars (already non-dimensionalized)
  rho0 = 1.0;
  c0 = 1.0;
  p0 = rho0*c0*c0/gma;
  u0 = xmach*cos(alpha/conv);
  v0 = xmach*sin(alpha/conv);
  ei0 = p0/((gma - 1.0)*rho0);
  et0 = rho0*(ei0 + 0.5*(u0*u0 + v0*v0));

  //now, read in the mesh file name before we init q, we have to know where to put it and bd cond
  char filename[bdim]; //to hold mesh file name
  printf("\nEnter name of mesh file with extension ->");
  fgets(buff,bdim,stdin); //read in mesh file name
  sscanf(buff,"%s",&filename); //place in variable

  //now, declare mesh object to hold unstructured mesh
  //all functions will be in its class definition so that vars needn't be passed in/out of functions 
  //this is okay since these are definitely global values
  Mesh mesh_obj;

  //place gma in mesh_obj var for use from here on out
  mesh_obj.gamma = gma;

  //now, call mesh_read to read in mesh
  mesh_obj.mesh_read(filename);

  //now, create hash table for compressed row storage routine
  mesh_obj.Lhash();
  
  //now, create compressed row storage arrays
  mesh_obj.crs();

  //now, read in all inviscid inner boundaries (walls)
  printf("\nIf there are no inviscid inner boundaries, then use -1 for both of the following indices!");
  printf("\nEnter inviscid wall inner boundary numbers (must have two), C++ indexed, with a space in between  ->");
  fgets(buff,bdim,stdin); //read in inner bd numbers
  sscanf(buff,"%d %d",&mesh_obj.ib1,&mesh_obj.ib2); //place in variables
  //error checking ..will be fine with -1 -1
  if (mesh_obj.ib1 >= mesh_obj.nb || mesh_obj.ib2 >= mesh_obj.nb)
    {
    printf("\nYou entered a boundary that doesn't exist.  Check your mesh file and be sure you used C++ indexing (decrement by 1).  Exiting....\n");
    fflush(stdout);
    exit(0);
    }

  //now, compute areas, compute len_scale, initialize q
  //pass in initial vals of prim vars
  //the passing of vars is myoot if using a restart file
  mesh_obj.mesh_init(rho0, u0, v0, et0);

  //finalize last few details before beginning linear solve
  //check for implicit or explicit solve
  printf("\nEnter 0 for explicit solve or 1 for implicit solve ->");
  fgets(buff,bdim,stdin); //read in impex
  sscanf(buff,"%d",&mesh_obj.impex); //place in variable
  //error checking
  if (mesh_obj.impex < 0 || mesh_obj.impex > 1)
    {
    printf("\nYou must choose 0 for explicit solve or 1 for implicit solve.  Exiting....\n");
    fflush(stdout);
    exit(0);
    }
  //read in CFL
  double CFLinit = 0.0, CFLmax = 0.0; //declare regardless of implicit or explicit, for compiler
  int ramp = 0; //declare regardless of implicit or explicit, for compiler
  if (mesh_obj.impex == 0)
    {
    printf("\nEnter desired CFL number ->");
    fgets(buff,bdim,stdin); //read in mesh_obj.CFL
    sscanf(buff,"%lf",&mesh_obj.CFL); //place in variable
    //error checking
    if (mesh_obj.CFL <= 0.0)
      {
      printf("\nCondition to satisfy: CFL > 0.0.  Exiting....\n");
      fflush(stdout);
      exit(0);
      }
    }
  if (mesh_obj.impex == 1)
    {
    //read in initial and max and set mesh_obj version in loop   
    printf("\nEnter desired initial CFL number ->");
    fgets(buff,bdim,stdin); //read in CFLinit
    sscanf(buff,"%lf",&CFLinit); //place in variable
    //error checking
    if (CFLinit <= 0.0)
      {
      printf("\nCondition to satisfy: CFL(initial) > 0.0.  Exiting....\n");
      fflush(stdout);
      exit(0);
      }
    printf("\nEnter desired maximum CFL number ->");
    fgets(buff,bdim,stdin); //read in CFLmax
    sscanf(buff,"%lf",&CFLmax); //place in variable
    //error checking
    if (CFLmax <= 0.0)
      {
      printf("\nCondition to satisfy: CFL(maximum) > 0.0.  Exiting....\n");
      fflush(stdout);
      exit(0);
      }
    //enter number of iterations to ramp up to max CFL
    printf("\nEnter number of iterations to ramp up to max CFL ->");
    fgets(buff,bdim,stdin); //read in ramp
    sscanf(buff,"%d",&ramp); //place in variable
    //error checking
    if (ramp <= 1)
      {
      printf("\nCondition to satisfy: ramp >= 2.  Exiting....\n");
      fflush(stdout);
      exit(0);
      }
    //read in order
    printf("\nEnter 1 for first order solve or 2 for second order solve ->");
    fgets(buff,bdim,stdin); //read in order
    sscanf(buff,"%d",&mesh_obj.order); //place in variable
    //error checking
    if (mesh_obj.order < 1 || mesh_obj.order > 2)
      {
      printf("\nMust choose 1 for first order or 2 for second order.  Exiting....\n");
      fflush(stdout);
      exit(0);
      }
    }

  //read in max number of iterations to perform before recompute (unless RMS low enough sooner)
  int itermax = 0;
  printf("\nEnter max number of iterations to be performed before recompute ->");
  fgets(buff,bdim,stdin); //read in itermax
  sscanf(buff,"%d",&itermax); //place in variable
  //error checking
  if (itermax <= 0)
    {
    printf("\nYou cannot solve a matrix in %d iterations!  Exiting....\n",itermax);
    fflush(stdout);
    exit(0);
    }

  //read in max number of iterations to perform of outer loop (unless RMS low enough sooner)
  int itertot = 0;
  printf("\nEnter max number of outer loop iterations to be performed ->");
  fgets(buff,bdim,stdin); //read in itermax
  sscanf(buff,"%d",&itertot); //place in variable
  //error checking
  if (itertot <= 0)
    {
    printf("\nYou cannot solve a matrix in %d iterations!  Exiting....\n",itertot);
    fflush(stdout);
    exit(0);
    }

  //read in omega (relaxation factor) for GS-SSOR
  printf("\nEnter relaxation factor ->");
  fgets(buff,bdim,stdin); //read in omega
  sscanf(buff,"%lf",&mesh_obj.omega); //place in variable
  //error checking
  if (mesh_obj.omega <= 0.0 || mesh_obj.omega > 1.0)
    {
    printf("\nCondition to satisfy: 0.0 < omega <=1.0.  Exiting....\n");
    fflush(stdout);
    exit(0);
    }

  //read in convergence value
  printf("\nEnter desired value for convergence ->");
  fgets(buff,bdim,stdin); //read in cvg
  sscanf(buff,"%lf",&mesh_obj.cvg); //place in variable
  //error checking
  if (mesh_obj.cvg > 1.0e-3 || mesh_obj.cvg < 1.0e-16)
    {
    printf("\nCondition to satisfy: 1.0e-16 <= cvg <=1.0e-3.  Exiting....\n");
    fflush(stdout);
    exit(0);
    }

  //init RMS error high for debugging
  double RMSerror = 1000.0;

  //init outer loop iteration counter
  int iter = 0;

  //init vars to hold iterations performed by inner loop and sum total of the same
  int tempiter = 0, tempitersum = 0;

  //set flag to 0, to allow our loop to continue until RHS matches
  int flag = 0;

  //declare ResidErr, init high so goes thru loop
  double ResidErr = 1000.0;

  //need to set pretty, typecast max CFL
  int printCFL = (int) CFLmax;

  // we print rms (resid err) to file with iterations as we loop through
  FILE *fp = NULL;

  buff[0]='\0';
  if (abs(alpha2) < 10 && alpha2 >= 0 && itermax < 100 && printCFL < 100)
    sprintf(buff,"RMSEulerIter_%d_+0%d_%1.2f_0%d_0%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 < 0 && itermax < 100 && printCFL < 100)
    sprintf(buff,"RMSEulerIter_%d_-0%d_%1.2f_0%d_0%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 >= 0 && printCFL < 100)
    sprintf(buff,"RMSEulerIter_%d_+0%d_%1.2f_%d_0%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 < 0 && printCFL < 100)
    sprintf(buff,"RMSEulerIter_%d_-0%d_%1.2f_%d_0%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 >= 0 && itermax < 100)
    sprintf(buff,"RMSEulerIter_%d_+0%d_%1.2f_0%d_%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 < 0 && itermax < 100)
    sprintf(buff,"RMSEulerIter_%d_-0%d_%1.2f_0%d_%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 >= 0)
    sprintf(buff,"RMSEulerIter_%d_+0%d_%1.2f_%d_%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 < 0)
    sprintf(buff,"RMSEulerIter_%d_-0%d_%1.2f_%d_%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (alpha2 < 0 && itermax < 100 && printCFL < 100)
    sprintf(buff,"RMSEulerIter_%d_%d_%1.2f_0%d_0%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 >= 0 && itermax < 100 && printCFL < 100)
    sprintf(buff,"RMSEulerIter_%d_+%d_%1.2f_0%d_0%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 >= 0 && itermax < 100)
    sprintf(buff,"RMSEulerIter_%d_+%d_%1.2f_0%d_%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 < 0 && itermax < 100)
    sprintf(buff,"RMSEulerIter_%d_%d_%1.2f_0%d_%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 >= 0 && printCFL < 100)
    sprintf(buff,"RMSEulerIter_%d_+%d_%1.2f_%d_0%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 < 0 && printCFL < 100)
    sprintf(buff,"RMSEulerIter_%d_%d_%1.2f_%d_0%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 < 0)
    sprintf(buff,"RMSEulerIter_%d_%d_%1.2f_%d_%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 >= 0)
    sprintf(buff,"RMSEulerIter_%d_+%d_%1.2f_%d_%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
    

  printf("\nImplicit Euler RMS Error (Iterations) Filename = <%s>\n",buff);
  
  // Open file for writing of rms errors
  if ((fp = fopen(buff,"w")) == 0)
    {
    printf("\nError opening file <%s>.",buff);
    exit(0);
    }

  // we print rms (resid err) to file with timing as we loop through
  /*FILE *fp2 = NULL;
 
  buff[0]='\0';
  strcat(buff,"RMSEulerTime.dat");
  printf("\nImplicit Euler RMS Error (Timing) Filename = <%s>\n",buff);
  
  // Open file for writing of rms errors
  if ((fp2 = fopen(buff,"w")) == 0)
    {
    printf("\nError opening file <%s>.",buff);
    exit(0);
    }*/ 

  //start timing..we will want time value from beginning at each interval
  gettimeofday(&start_time,NULL);

  //now, we begin our outer linear solve loop
  do
    {
    iter++; //increment outer loop iteration counter

    if (mesh_obj.impex == 1 && mesh_obj.order == 2)
      {
      //now, find dq/dx and dq/dy
      mesh_obj.gradients();
      }

    if (mesh_obj.impex == 1)
      {
      //ramp up CFL number
      //cast ramp and iter as doubles for calc
      mesh_obj.CFL = CFLinit + ((((double)(iter)) - 1.0)/(((double)(ramp)) - 1.0)) * (CFLmax - CFLinit);

      //if CFL goes over max, then reset to max
      if (mesh_obj.CFL > CFLmax)
        mesh_obj.CFL = CFLmax;
 
      //printf("\nCFL this run  = %lf\n",mesh_obj.CFL);
      //wait(3); //debug
      }

    //now, create A matrix (mat)
    mesh_obj.mat_build();

    //now, compute residual (RHS) and fill in A matrix with dR/dQ
    mesh_obj.residual();

    //now, check and see if RHS is zero vector
    ResidErr = 0.0; //reset to begin summing

    for (i = 0; i < mesh_obj.nn; i++)
      for (j = 0; j < 4; j++)
        ResidErr += (mesh_obj.RHS[i][j])*(mesh_obj.RHS[i][j]);

    //finalize error
    ResidErr = sqrt(ResidErr/(4*mesh_obj.nn));

    //if it has not changed, set flag to skip next solve, since we will get no change
    if (ResidErr <= 1.0e-15)
      flag = 1; 

    //only run through if RHS has changed (flag == 0)
    if (!flag)
      {
      //set tempiter to itermax to be passed in and then out by ref as iter performed
      tempiter = itermax;

      //now, do linear solve (will pass out iterations, and return RMSerror)
      RMSerror = mesh_obj.lin_solve(tempiter);

      //begin summing total internal iterations performed
      tempitersum += tempiter;

      //print current iter and RMS for itermediate steps only
      //printf("\nCurrent RMSerror between R(n) and R(n+1) of %16.10e after %d external iterations and a sum of %d internal iterations.\n",ResidErr,iter,tempitersum);
      //wait(3); //debug
      }

    //print rms value and iteration number to one file
    fprintf(fp,"%d %19.10e\n",iter,ResidErr);

    //stop timing
    gettimeofday(&end_time,NULL);
    delta_time = (double)(end_time.tv_sec - start_time.tv_sec);
    delta_time += (double)(end_time.tv_usec - start_time.tv_usec)/1.0e6;

    //print rms value and time to different file
    //fprintf(fp2,"%19.10e %19.10e\n",delta_time,ResidErr);

    } while (ResidErr > 1.0e-15 && iter < itertot);

  fclose(fp); //close file
  //fclose(fp2); //close file2
  fp = NULL;  //reset pointer
  //fp2 = NULL;  //reset pointer2

  //print final RMS and iteration info
  printf("\nYou converged to an RMSerror between R(n) and R(n+1) of %16.10e in %d external iterations and %d total internal iterations.\n",ResidErr,iter,tempitersum);

  //only print cp on inner bd, do separate files..XMGR style
  for (k = 0; k < 2; k++)
    {
    switch(k)
      {
      case 0:
        l = mesh_obj.ib1;    
      break;
      case 1:
        l = mesh_obj.ib2;     
      break;
      default:
        printf("\nYou are trying to print cp on more than two boundaries!\n");
        fflush(stdout);
        exit(0);    
      break;
      }
    //calc cp and print to file
    buff[0]='\0';
  if (abs(alpha2) < 10 && alpha2 >= 0 && itermax < 100 && printCFL < 100)
    sprintf(buff,"CPEuler_%d_%d_+0%d_%1.2f_0%d_0%d.dat",l,mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 < 0 && itermax < 100 && printCFL < 100)
    sprintf(buff,"CPEuler_%d_%d_-0%d_%1.2f_0%d_0%d.dat",l,mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 >= 0 && printCFL < 100)
    sprintf(buff,"CPEuler_%d_%d_+0%d_%1.2f_%d_0%d.dat",l,mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 < 0 && printCFL < 100)
    sprintf(buff,"CPEuler_%d_%d_-0%d_%1.2f_%d_0%d.dat",l,mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 >= 0 && itermax < 100)
    sprintf(buff,"CPEuler_%d_%d_+0%d_%1.2f_0%d_%d.dat",l,mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 < 0 && itermax < 100)
    sprintf(buff,"CPEuler_%d_%d_-0%d_%1.2f_0%d_%d.dat",l,mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 >= 0)
    sprintf(buff,"CPEuler_%d_%d_+0%d_%1.2f_%d_%d.dat",l,mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 < 0)
    sprintf(buff,"CPEuler_%d_%d_-0%d_%1.2f_%d_%d.dat",l,mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (alpha2 < 0 && itermax < 100 && printCFL < 100)
    sprintf(buff,"CPEuler_%d_%d_%d_%1.2f_0%d_0%d.dat",l,mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 >= 0 && itermax < 100 && printCFL < 100)
    sprintf(buff,"CPEuler_%d_%d_+%d_%1.2f_0%d_0%d.dat",l,mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 >= 0 && itermax < 100)
    sprintf(buff,"CPEuler_%d_%d_+%d_%1.2f_0%d_%d.dat",l,mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 < 0 && itermax < 100)
    sprintf(buff,"CPEuler_%d_%d_%d_%1.2f_0%d_%d.dat",l,mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 >= 0 && printCFL < 100)
    sprintf(buff,"CPEuler_%d_%d_+%d_%1.2f_%d_0%d.dat",l,mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 < 0 && printCFL < 100)
    sprintf(buff,"CPEuler_%d_%d_%d_%1.2f_%d_0%d.dat",l,mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 < 0)
    sprintf(buff,"CPEuler_%d_%d_%d_%1.2f_%d_%d.dat",l,mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 >= 0)
    sprintf(buff,"CPEuler_%d_%d_+%d_%1.2f_%d_%d.dat",l,mesh_obj.order,alpha2,xmach,itermax,printCFL);

    printf("\nImplicit Euler Pressure Coefficient XMGR Filename = <%s>\n",buff);
  
    // Open file for writing of cp
    if ((fp = fopen(buff,"w")) == 0)
      { 
      printf("\nError opening file <%s>.",buff);
      exit(0);
      }

    //we will use i and j to be sure we only print out each bd node once
    for (i = 0; i < mesh_obj.nbs[l]; i++)
      {
      for (j = 0; j < 2; j++)
        {
        if (i == 0 || j == 1)
          {
          //use dummy index to keep cleaner
          n = mesh_obj.bs[l][i][j];
          //use dummy var to keep cleaner
          p = (gma - 1.0)*(mesh_obj.q[n][3] - 0.5*((mesh_obj.q[n][1]*mesh_obj.q[n][1] + mesh_obj.q[n][2]*mesh_obj.q[n][2])/mesh_obj.q[n][0]));
          cp = (2.0*(gma*p - 1.0))/(gma*xmach*xmach);
          x = mesh_obj.nodep[n][0];
          y = mesh_obj.nodep[n][1];
          //print to file x, cp, y
          fprintf(fp,"%19.10e %19.10e %19.10e\n",x,cp,y);
          }
        }
      }
    fclose(fp);  //close file
    fp = NULL;  //reset pointer
    }

  //tecplot file of u, v, rho, et, p, cp  
  buff[0]='\0';

  if (abs(alpha2) < 10 && alpha2 >= 0 && itermax < 100 && printCFL < 100)
    sprintf(buff,"TecplotEuler_%d_+0%d_%1.2f_0%d_0%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 < 0 && itermax < 100 && printCFL < 100)
    sprintf(buff,"TecplotEuler_%d_-0%d_%1.2f_0%d_0%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 >= 0 && printCFL < 100)
    sprintf(buff,"TecplotEuler_%d_+0%d_%1.2f_%d_0%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 < 0 && printCFL < 100)
    sprintf(buff,"TecplotEuler_%d_-0%d_%1.2f_%d_0%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 >= 0 && itermax < 100)
    sprintf(buff,"TecplotEuler_%d_+0%d_%1.2f_0%d_%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 < 0 && itermax < 100)
    sprintf(buff,"TecplotEuler_%d_-0%d_%1.2f_0%d_%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 >= 0)
    sprintf(buff,"TecplotEuler_%d_+0%d_%1.2f_%d_%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (abs(alpha2) < 10 && alpha2 < 0)
    sprintf(buff,"TecplotEuler_%d_-0%d_%1.2f_%d_%d.dat",mesh_obj.order,abs(alpha2),xmach,itermax,printCFL);
  else if (alpha2 < 0 && itermax < 100 && printCFL < 100)
    sprintf(buff,"TecplotEuler_%d_%d_%1.2f_0%d_0%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 >= 0 && itermax < 100 && printCFL < 100)
    sprintf(buff,"TecplotEuler_%d_+%d_%1.2f_0%d_0%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 >= 0 && itermax < 100)
    sprintf(buff,"TecplotEuler_%d_+%d_%1.2f_0%d_%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 < 0 && itermax < 100)
    sprintf(buff,"TecplotEuler_%d_%d_%1.2f_0%d_%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 >= 0 && printCFL < 100)
    sprintf(buff,"TecplotEuler_%d_+%d_%1.2f_%d_0%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 < 0 && printCFL < 100)
    sprintf(buff,"TecplotEuler_%d_%d_%1.2f_%d_0%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 < 0)
    sprintf(buff,"TecplotEuler_%d_%d_%1.2f_%d_%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);
  else if (alpha2 >= 0)
    sprintf(buff,"TecplotEuler_%d_+%d_%1.2f_%d_%d.dat",mesh_obj.order,alpha2,xmach,itermax,printCFL);

  printf("\nImplicit Euler Tecplot Filename = <%s>\n",buff);
  
  // Open file for writing
  if ((fp = fopen(buff,"w")) == 0)
    { 
    printf("\nError opening file <%s>.",buff);
    exit(0);
    }

  //print techplot vars
  fprintf(fp, "TITLE=\"Implicit Euler\"\n");
  fprintf(fp,"VARIABLES=\"X\", \"Y\", \"Cp\", \"rho\", \"u\", \"v\", \"p\", \"et\"\n");
  fprintf(fp,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n", mesh_obj.nn, mesh_obj.nt);

  //print to file
  for (i = 0; i < mesh_obj.nn; i++)
    {
    //use dummy var to keep cleaner
    p = (gma - 1.0)*(mesh_obj.q[i][3] - 0.5*((mesh_obj.q[i][1]*mesh_obj.q[i][1] + mesh_obj.q[i][2]*mesh_obj.q[i][2])/mesh_obj.q[i][0]));
    cp = (2.0*(gma*p - 1.0))/(gma*xmach*xmach);
    x = mesh_obj.nodep[i][0];
    y = mesh_obj.nodep[i][1];
    fprintf(fp, "%16.10lf %16.10lf %16.10lf %16.10lf %16.10lf %16.10lf %16.10lf %16.10lf\n", x, y, cp, mesh_obj.q[i][0], mesh_obj.q[i][1]/mesh_obj.q[i][0], mesh_obj.q[i][2]/mesh_obj.q[i][0], p, mesh_obj.q[i][3]);
    }

  fprintf(fp,"\n");

  for (i = 0; i < mesh_obj.nt; i++)
    {
    //use dummy var to keep cleaner
    n0 = mesh_obj.tri[i][0];
    n1 = mesh_obj.tri[i][1];
    n2 = mesh_obj.tri[i][2];
    fprintf(fp, "%d %d %d \n", n0 + 1, n1 + 1, n2 + 1);
    }      

  fclose(fp);

  //need to print out restart mesh file to be read in if we desire to start first order and then go second
  mesh_obj.mesh_out(filename);

  //debug
  //printf("\nMust be a destructor problem.");
  //fflush(stdout);
  
  //mesh_obj not a pointer, so it auto-destructs

  //check for nan, so go back, run first order, restart
  if (isnan(ResidErr))
    return(0);
  else if (ResidErr > 1.0e-12)
    return(-1);
  else
    return(1);
  }
