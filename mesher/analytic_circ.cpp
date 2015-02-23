#include "analytic_circ.h"

//takes coords and returns function value at node

double analytic_circ(double x, double y)
  {
  //first, perturb point
  x -= 0.5;
  y -= 0.5;

  //now, determine length between it and 0,0
  double dr = sqrt(x*x + y*y);

  //check and see if inside ring else, outside ring
  if (dr <= 0.25 && dr >= 0.125)
    return (1.0);
  else
    return (0.0);
  } 
