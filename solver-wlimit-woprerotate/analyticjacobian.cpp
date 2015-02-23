#include "analyticjacobian.h"

//NOTE: The special linearization for inner, inviscid boundary is in residual.cpp! 

void analyticjacobian(double *ql, double *qr, double gamma, double xnorm, double ynorm, double rlen, double *fp, double *fm, double **afp, double **afm)
  {
  //declare variables
  double gm1, gp1, gm1g, gp1g, ggm1; //all gamma's common variations
  int i, j; //counters
  double rhol, ul, vl, pl, ubarl, cl, etl, rhor, ur, vr, pr, ubarr, cr, etr; //primitive vars to be differentiated w.r.t. qr and ql
  //primitive vars after differentiation
  double rhol_q1, ul_q1, vl_q1, pl_q1, ubarl_q1, cl_q1, rhol_q2, ul_q2, vl_q2, pl_q2, ubarl_q2, cl_q2, rhol_q3, ul_q3, vl_q3, pl_q3, ubarl_q3, cl_q3, rhol_q4, ul_q4, vl_q4, pl_q4, ubarl_q4, cl_q4, rhor_q1, ur_q1, vr_q1, pr_q1, ubarr_q1, cr_q1, rhor_q2, ur_q2, vr_q2, pr_q2, ubarr_q2, cr_q2, rhor_q3, ur_q3, vr_q3, pr_q3, ubarr_q3, cr_q3, rhor_q4, ur_q4, vr_q4, pr_q4, ubarr_q4, cr_q4; 

  //initialize gamma parameters
  gm1 = gamma - 1.0;
  gp1 = gamma + 1.0;
  gm1g = gm1/gamma;
  gp1g = gp1/gamma;
  ggm1 = gamma*gm1;

  //ANALYTIC JACOBIANS
  //now, retrieve primitive vars from ql and qr
  rhol = ql[0];
  ul = ql[1]/ql[0];
  vl = ql[2]/ql[0];
  pl = gm1*(ql[3] - 0.5*((ql[1]*ql[1] + ql[2]*ql[2])/(ql[0])));
  ubarl = xnorm*ul + ynorm*vl;
  cl = sqrt((gamma*pl)/rhol);
  etl = ql[3];

  rhor = qr[0];
  ur = qr[1]/qr[0];
  vr = qr[2]/qr[0];
  pr = gm1*(qr[3] - 0.5*((qr[1]*qr[1] + qr[2]*qr[2])/(qr[0])));
  ubarr = xnorm*ur + ynorm*vr;
  cr = sqrt((gamma*pr)/rhor);
  etr = qr[3];

  //create arrays of primitive vars w.r.t. each element of qr and ql (for easy creation of flux jacobian)
  double rhol_ql[4], ul_ql[4], vl_ql[4], pl_ql[4], ubarl_ql[4], cl_ql[4], el_ql[4], rhor_qr[4], ur_qr[4], vr_qr[4], pr_qr[4], ubarr_qr[4], cr_qr[4], er_qr[4];
  
  //fill-in values of partials (using the definition of variable in qr and ql to find partial)
  rhol_ql[0] = 1.0;
  rhol_ql[1] = 0.0;
  rhol_ql[2] = 0.0;
  rhol_ql[3] = 0.0;

  ul_ql[0] = -ql[1]/(ql[0]*ql[0]);
  ul_ql[1] = 1.0/ql[0];
  ul_ql[2] = 0.0;
  ul_ql[3] = 0.0;

  vl_ql[0] = -ql[2]/(ql[0]*ql[0]);
  vl_ql[1] = 0.0;
  vl_ql[2] = 1.0/ql[0];
  vl_ql[3] = 0.0;

  pl_ql[0] = (0.5*gm1*(ql[1]*ql[1] + ql[2]*ql[2]))/(ql[0]*ql[0]);
  pl_ql[1] = (1.0 - gamma)*ql[1]/ql[0];
  pl_ql[2] = (1.0 - gamma)*ql[2]/ql[0];
  pl_ql[3] = gm1;

  ubarl_ql[0] = (-xnorm*ql[1] - ynorm*ql[2])/(ql[0]*ql[0]);
  ubarl_ql[1] = xnorm/ql[0];
  ubarl_ql[2] = ynorm/ql[0];
  ubarl_ql[3] = 0.0;

  el_ql[0] = 0.0;
  el_ql[1] = 0.0;
  el_ql[2] = 0.0;
  el_ql[3] = 1.0;

  cl_ql[0] = (1.0/(2.0*sqrt(ggm1*ql[3]/ql[0] - 0.5*ggm1*(ql[1]*ql[1] + ql[2]*ql[2])/(ql[0]*ql[0]))))*((-ggm1*ql[3])/(ql[0]*ql[0]) + ((0.5*ggm1)*(ql[1]*ql[1] + ql[2]*ql[2])*2.0*ql[0])/(ql[0]*ql[0]*ql[0]*ql[0]));
  cl_ql[1] = (1.0/(2.0*sqrt(ggm1*ql[3]/ql[0] - 0.5*ggm1*(ql[1]*ql[1] + ql[2]*ql[2])/(ql[0]*ql[0]))))*((-(ggm1)/(ql[0]*ql[0]))*ql[1]);
  cl_ql[2] = (1.0/(2.0*sqrt(ggm1*ql[3]/ql[0] - 0.5*ggm1*(ql[1]*ql[1] + ql[2]*ql[2])/(ql[0]*ql[0]))))*((-(ggm1)/(ql[0]*ql[0]))*ql[2]);
  cl_ql[3] = (1.0/(2.0*sqrt(ggm1*ql[3]/ql[0] - 0.5*ggm1*(ql[1]*ql[1] + ql[2]*ql[2])/(ql[0]*ql[0]))))*((ggm1)/ql[0]);

  rhor_qr[0] = 1.0;
  rhor_qr[1] = 0.0;
  rhor_qr[2] = 0.0;
  rhor_qr[3] = 0.0;

  ur_qr[0] = -qr[1]/(qr[0]*qr[0]);
  ur_qr[1] = 1.0/qr[0];
  ur_qr[2] = 0.0;
  ur_qr[3] = 0.0;

  vr_qr[0] = -qr[2]/(qr[0]*qr[0]);
  vr_qr[1] = 0.0;
  vr_qr[2] = 1.0/qr[0];
  vr_qr[3] = 0.0;

  pr_qr[0] = (0.5*gm1*(qr[1]*qr[1] + qr[2]*qr[2]))/(qr[0]*qr[0]);
  pr_qr[1] = (1.0 - gamma)*qr[1]/qr[0];
  pr_qr[2] = (1.0 - gamma)*qr[2]/qr[0];
  pr_qr[3] = gm1;

  ubarr_qr[0] = (-xnorm*qr[1] - ynorm*qr[2])/(qr[0]*qr[0]);
  ubarr_qr[1] = xnorm/qr[0];
  ubarr_qr[2] = ynorm/qr[0];
  ubarr_qr[3] = 0.0;

  er_qr[0] = 0.0;
  er_qr[1] = 0.0;
  er_qr[2] = 0.0;
  er_qr[3] = 1.0;

  cr_qr[0] = (1.0/(2.0*sqrt(ggm1*qr[3]/qr[0] - 0.5*ggm1*(qr[1]*qr[1] + qr[2]*qr[2])/(qr[0]*qr[0]))))*((-ggm1*qr[3])/(qr[0]*qr[0]) + ((0.5*ggm1)*(qr[1]*qr[1] + qr[2]*qr[2])*2.0*qr[0])/(qr[0]*qr[0]*qr[0]*qr[0]));
  cr_qr[1] = (1.0/(2.0*sqrt(ggm1*qr[3]/qr[0] - 0.5*ggm1*(qr[1]*qr[1] + qr[2]*qr[2])/(qr[0]*qr[0]))))*((-(ggm1)/(qr[0]*qr[0]))*qr[1]);
  cr_qr[2] = (1.0/(2.0*sqrt(ggm1*qr[3]/qr[0] - 0.5*ggm1*(qr[1]*qr[1] + qr[2]*qr[2])/(qr[0]*qr[0]))))*((-(ggm1)/(qr[0]*qr[0]))*qr[2]);
  cr_qr[3] = (1.0/(2.0*sqrt(ggm1*qr[3]/qr[0] - 0.5*ggm1*(qr[1]*qr[1] + qr[2]*qr[2])/(qr[0]*qr[0]))))*((ggm1)/qr[0]);
  

  int flag1p = 0, flag2p = 0, flag1m = 0, flag2m = 0; //check for bidirectional (incorrect) flow
  
  //adjust for cases with supersonic flow
  if (ubarl/cl > 1.0)
    {
    //now, fill in afp
    for (i = 0; i < 4; i++)
      {
      afp[0][i] = rlen*(rhol_ql[i]*ubarl + ubarl_ql[i]*rhol);
      afp[1][i] = rlen*(rhol_ql[i]*ul*ubarl + ul_ql[i]*rhol*ubarl + ubarl_ql[i]*rhol*ul + xnorm*pl_ql[i]); 
      afp[2][i] = rlen*(rhol_ql[i]*vl*ubarl + vl_ql[i]*rhol*ubarl + ubarl_ql[i]*rhol*vl + ynorm*pl_ql[i]);
      afp[3][i] = rlen*(el_ql[i]*ubarl + ubarl_ql[i]*etl + pl_ql[i]*ubarl + ubarl_ql[i]*pl);
      }

    flag1p == 1;
    }

  else if (ubarl/cl < -1.0)
    {
    //now, fill in afp
    for (i = 0; i < 4; i++)
      {
      afp[0][i] = 0.0;
      afp[1][i] = 0.0; 
      afp[2][i] = 0.0;
      afp[3][i] = 0.0;
      }

    flag2p == 1;
    }
 
  else
    {
    //now, fill in afp
    for (i = 0; i < 4; i++)
      {
      afp[0][i] = 0.25*(rhol_ql[i]*cl*(ubarl/cl + 1.0)*(ubarl/cl + 1.0) + cl_ql[i]*rhol*(ubarl/cl + 1.0)*(ubarl/cl + 1.0) + 2.0*rhol*cl*(ubarl/cl + 1.0)*(1.0/(cl*cl))*(ubarl_ql[i]*cl - cl_ql[i]*ubarl))*rlen;
      afp[1][i] = fp[0]*((xnorm/gamma)*(-ubarl_ql[i] + 2.0*cl_ql[i]) + ul_ql[i]) + ((xnorm/gamma)*(-ubarl + 2.0*cl) + ul)*afp[0][i]; 
      afp[2][i] = fp[0]*((ynorm/gamma)*(-ubarl_ql[i] + 2.0*cl_ql[i]) + vl_ql[i]) + ((ynorm/gamma)*(-ubarl + 2.0*cl) + vl)*afp[0][i];
      afp[3][i] = fp[0]*((-2.0/gp1)*ubarl*ubarl_ql[i] + (2.0/gp1)*(ubarl*cl_ql[i] + cl*ubarl_ql[i]) + (4.0/(gamma*gamma - 1.0))*cl*cl_ql[i] + ul*ul_ql[i] + vl*vl_ql[i]) + ((-gm1*ubarl*ubarl + 2.0*gm1*cl*ubarl + 2.0*cl*cl)/(gamma*gamma - 1.0) + 0.5*(ul*ul + vl*vl))*afp[0][i];
      }
    }

  if (ubarr/cr < -1.0)
    {
    //now, fill in afm
    for (i = 0; i < 4; i++)
      {
      afm[0][i] = rlen*(rhor_qr[i]*ubarr + ubarr_qr[i]*rhor);
      afm[1][i] = rlen*(rhor_qr[i]*ur*ubarr + ur_qr[i]*rhor*ubarr + ubarr_qr[i]*rhor*ur + xnorm*pr_qr[i]); 
      afm[2][i] = rlen*(rhor_qr[i]*vr*ubarr + vr_qr[i]*rhor*ubarr + ubarr_qr[i]*rhor*vr + ynorm*pr_qr[i]);
      afm[3][i] = rlen*(er_qr[i]*ubarr + ubarr_qr[i]*etr + pr_qr[i]*ubarr + ubarr_qr[i]*pr);
      }

    flag2m == 1;
    }

  else if (ubarr/cr > 1.0)
    {
    //now, fill in afm
    for (i = 0; i < 4; i++)
      {
      afm[0][i] = 0.0;
      afm[1][i] = 0.0; 
      afm[2][i] = 0.0;
      afm[3][i] = 0.0;
      }

    flag1m == 1;
    }

  else 
    {
    //now, fill in afm
    for (i = 0; i < 4; i++)
      {
      afm[0][i] = -0.25*(rhor_qr[i]*cr*(ubarr/cr - 1.0)*(ubarr/cr - 1.0) + cr_qr[i]*rhor*(ubarr/cr - 1.0)*(ubarr/cr - 1.0) + 2.0*rhor*cr*(ubarr/cr - 1.0)*(1.0/(cr*cr))*(ubarr_qr[i]*cr - cr_qr[i]*ubarr))*rlen;
      afm[1][i] = fm[0]*((xnorm/gamma)*(-ubarr_qr[i] - 2.0*cr_qr[i]) + ur_qr[i]) + ((xnorm/gamma)*(-ubarr - 2.0*cr) + ur)*afm[0][i]; 
      afm[2][i] = fm[0]*((ynorm/gamma)*(-ubarr_qr[i] - 2.0*cr_qr[i]) + vr_qr[i]) + ((ynorm/gamma)*(-ubarr - 2.0*cr) + vr)*afm[0][i];
      afm[3][i] = fm[0]*((-2.0/gp1)*ubarr*ubarr_qr[i] - (2.0/gp1)*(ubarr*cr_qr[i] + cr*ubarr_qr[i]) + (4.0/(gamma*gamma - 1.0))*cr*cr_qr[i] + ur*ur_qr[i] + vr*vr_qr[i]) + ((-gm1*ubarr*ubarr - 2.0*gm1*cr*ubarr + 2.0*cr*cr)/(gamma*gamma - 1.0) + 0.5*(ur*ur + vr*vr))*afm[0][i];
      }
    }

  //debug, check for impossible flow
  if (flag1m == 1 && flag1p == 1)
    {
    printf("\nYour flow has a split personality.");
    fflush(stdout);
    exit(0);
    }

  if (flag2m == 1 && flag2p == 1)
    {
    printf("\nYour flow has a split personality.");
    fflush(stdout);
    exit(0);
    }

  return;
  }
