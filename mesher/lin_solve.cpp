#include "lin_solve.h"

//lin_solve takes the number of internal loop iterations, number of nodes to loop thru, mat (rebuilt in mat_build), tag to avoid solving at boundries, u and v which are either physical nodes (winslow) or changes in location (L-E), convergence value to jump out at, relaxation factor, and compressed row storage index arrays
//it will pass out RMSerror for convergence check
double lin_solve(int &iterin, int nn, double ***mat, int *tag, double *u, double *v, double cvg, double omega, int *ia, int *ja, int *iau)
  {
  int i, j, flag; //counters and flag
  double delu, delv, ustar, vstar; //delta u and delta v, placeholders for new u and v

  double *RHS; //initialize temporary RHS side array
  RHS = (double*)calloc(2,sizeof(double));

  //set RMSerr high to assure we start through
  double RMSerr = 1000.0;

  //start iter at 0 each time through
  int iter = 0;

  //now, we perform our solve...we do this a prescribed number or times, then we will go back, update, and do again
  //however, if RMSerr gets low enough we will jump out of this and it will kick out of main loop too!
  do
    {
    iter++; //increment iteration counter
    RMSerr = 0.0; //reset RMSerr to begin summing
    flag = 0; //reset flag keeping up with bd nodes

    //now, we loop through nodes and do solve...we are updating u and v as we go, so this is gauss-seidel
    //in order to rectify bias situations, we will do the loop in the opposite direction every other time thru  
    if (iter % 2 != 0)
      {
      for (i = 0; i < nn; i++)
        {
        if (tag[i] == 0)
          {
          flag++; //count up number of bd points to take away from nn when computing RMS
          continue; //we don't do this for boundary nodes!
          }

        //init RHS to 0
        for (j = 0; j < 2; j++)
          RHS[j] = 0.0;

        //now, transfer off-diag elements times u and v to RHS
        for (j = ia[i]; j < ia[i+1]; j++)
          {
          if (ja[j] == i)
            continue; //this is the diagonal element
          //now, transfer each over changing signs and "summing"
          //again, even though the off diag are 0 in winslow, we make things general
          RHS[0] -= mat[j][0][0]*u[ja[j]] + mat[j][0][1]*v[ja[j]];
          RHS[1] -= mat[j][1][1]*v[ja[j]] + mat[j][1][0]*u[ja[j]];
          }

        //use ustar and vstar as holders for the updated values of u and v before replacing them so we can check delu and delv and RMS
        ustar = RHS[0]*mat[iau[i]][0][0] + RHS[1]*mat[iau[i]][0][1];
        vstar = RHS[0]*mat[iau[i]][1][0] + RHS[1]*mat[iau[i]][1][1];

        //now, obtain delu and delv
        delu = ustar - u[i];
        delv = vstar - v[i];

        //now, replace u and v with relaxation to make G-S with SSOR
        u[i] += omega*delu;
        v[i] += omega*delv;
    
        //start summing RMSerr
        RMSerr += delu*delu+delv*delv;
        }
      }

    if (iter % 2 == 0)
      {
      for (i = nn-1; i > -1; i--)
        {
        if (tag[i] == 0)
          {
          flag++; //count up number of bd points to take away from nn when computing RMS
          continue; //we don't do this for boundary nodes!
          }

        //init RHS to 0
        for (j = 0; j < 2; j++)
          RHS[j] = 0.0;

        //now, transfer off-diag elements times u and v to RHS
        for (j = ia[i]; j < ia[i+1]; j++)
          {
          if (ja[j] == i)
            continue; //this is the diagonal element
          //now, transfer each over changing signs and "summing"
          //again, even though the off diag are 0 in winslow, we make things general
          RHS[0] -= mat[j][0][0]*u[ja[j]] + mat[j][0][1]*v[ja[j]];
          RHS[1] -= mat[j][1][1]*v[ja[j]] + mat[j][1][0]*u[ja[j]];
          }

        //use ustar and vstar as holders for the updated values of u and v before replacing them so we can check delu and delv and RMS
        ustar = RHS[0]*mat[iau[i]][0][0] + RHS[1]*mat[iau[i]][0][1];
        vstar = RHS[0]*mat[iau[i]][1][0] + RHS[1]*mat[iau[i]][1][1];

        //now, obtain delu and delv
        delu = ustar - u[i];
        delv = vstar - v[i];

        //now, replace u and v with relaxation to make G-S with SSOR
        u[i] += omega*delu;
        v[i] += omega*delv;
   
        //start summing RMSerr
        RMSerr += delu*delu+delv*delv;
        }
      }
    
    //finally, compute RMSerr
    RMSerr = sqrt(RMSerr/(nn-flag));

    //printf("\nYour RMS error is %16.10e\n",RMSerr);
   
    }while (iter < iterin && RMSerr > cvg);

  //printf("\nYou performed %d internal iterations\n",iter);

  //set iterin as iter for passing out in linear_elastic loop (will make accomodations in winslow for pasing by reference not causing detrimental overwrite)
  iterin = iter;

  freenull(RHS); //clean up memory allocated in the routine

  return(RMSerr);
  } 
