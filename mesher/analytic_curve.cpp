#include "analytic_curve.h"

//takes coords and returns function value at node

double analytic_curve(double x, double y)
  {
  //could use M_PI also
  double pi = 4.0*atan(1.0);

  //declare angle and cosine relational pos vars
  double ang, c;

  //now, compute angle (rad)
  ang = x*pi*2.0;

  //now, compute cos relational pos
  c = cos(ang)*0.25 + 0.5;

  //check and see if under curve, else above it
  if (y <= c)
    return (0.0);
  else
    return (1.0);
  } 
