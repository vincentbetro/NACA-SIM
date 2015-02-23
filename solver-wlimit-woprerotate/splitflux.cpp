#include "splitflux.h"

void splitflux(double *ql, double *qr, double gamma, double xnorm, double ynorm, double rlen, double *fps, double *fms)
  {
  //now, retrieve primitive vars from ql
  double rhol = ql[0];
  double ul = ql[1]/ql[0];
  double vl = ql[2]/ql[0];
  double pl = (gamma - 1.0)*(ql[3] - 0.5*((ql[1]*ql[1] + ql[2]*ql[2])/(ql[0])));
  double ubarl = xnorm*ul + ynorm*vl;
  double cl = sqrt((gamma*pl)/rhol);

  //now, retrieve primitive vars from qr
  double rhor = qr[0];
  double ur = qr[1]/qr[0];
  double vr = qr[2]/qr[0];
  double pr = (gamma - 1.0)*(qr[3] - 0.5*((qr[1]*qr[1] + qr[2]*qr[2])/(qr[0])));
  double ubarr = xnorm*ur + ynorm*vr;
  double cr = sqrt((gamma*pr)/rhor);

  int i;

  int flag1p = 0, flag2p = 0, flag1m = 0, flag2m = 0; //check for bidirectional (incorrect) flow
  
  //determine F+,F-, since for ubar/c >= 1.0, F+=F and for ubar/c <= -1.0, F-=F
  //use averages of right and left states
  if (ubarl/cl > 1.0)
    {
    fps[0] = rhol*ubarl; 
    fps[1] = rhol*ul*ubarl + xnorm*pl; 
    fps[2] = rhol*vl*ubarl + ynorm*pl;
    fps[3] = ubarl*(ql[3] + pl);

    flag1p == 1;
    }

  else if (ubarl/cl < -1.0)
    {
    fps[0] = 0.0; 
    fps[1] = 0.0; 
    fps[2] = 0.0;
    fps[3] = 0.0;

    flag2p == 1;
    }
 
  else
    {
    fps[0] = 0.25*rhol*cl*(ubarl/cl + 1.0)*(ubarl/cl + 1.0); 
    fps[1] = fps[0]*((xnorm/gamma)*(-ubarl + 2.0*cl) + ul); 
    fps[2] = fps[0]*((ynorm/gamma)*(-ubarl + 2.0*cl) + vl);
    fps[3] = fps[0]*((-(gamma - 1.0)*ubarl*ubarl + 2.0*(gamma - 1.0)*cl*ubarl + 2.0*cl*cl)/(gamma*gamma - 1.0) + 0.5*(ul*ul + vl*vl));
    }

  if (ubarr/cr < -1.0)
    {
    fms[0] = rhor*ubarr; 
    fms[1] = rhor*ur*ubarr + xnorm*pr; 
    fms[2] = rhor*vr*ubarr + ynorm*pr;
    fms[3] = ubarr*(qr[3] + pr);

    flag2m == 1;
    }

  else if (ubarr/cr > 1.0)
    {
    fms[0] = 0.0; 
    fms[1] = 0.0; 
    fms[2] = 0.0;
    fms[3] = 0.0;

    flag1m == 1;
    }

  else 
    {
    fms[0] = -0.25*rhor*cr*(ubarr/cr - 1.0)*(ubarr/cr - 1.0); 
    fms[1] = fms[0]*((xnorm/gamma)*(-ubarr - 2.0*cr) + ur); 
    fms[2] = fms[0]*((ynorm/gamma)*(-ubarr - 2.0*cr) + vr);
    fms[3] = fms[0]*((-(gamma - 1.0)*ubarr*ubarr - 2.0*(gamma - 1.0)*cr*ubarr + 2.0*cr*cr)/(gamma*gamma - 1.0) + 0.5*(ur*ur + vr*vr));
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

  //multiply F+, F- by length
  for (i = 0; i < 4; i++)
    {
    fps[i] *= rlen;
    fms[i] *= rlen;
    }

  return;
  }
