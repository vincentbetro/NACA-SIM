#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "refine_tri.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

double distance( double p1[2], double p2[2])
{
  double dx = p2[0]-p1[0];
  double dy = p2[1]-p1[1];
  double mag = sqrt(dx*dx+dy*dy);
  return mag;
}

double triangle_aspect_ratio(double p1[2], double p2[2], double p3[2])
{
  double a, b, c, s;
  double asp;

  a = distance(p1,p2);
  b = distance(p2,p3);
  c = distance(p3,p1);

  s = 0.5*(a+b+c);
  asp = (a*b*c)/(8.0*(s-a)*(s-b)*(s-c));

  return(asp);
}

int refine_tri(int conn[6], int cdim, int tri[][3], double **nodep)
{
  double asp1, asp2;

  // This routine returns the number of sub-triangles (at least 1)
  // The sub-triangle connectivity will be stored in the "tri" array

  assert(cdim >= 4);

  int nt = 0;

  int mask;
  mask = 0;
  if (conn[3] >= 0) mask += 1;
  if (conn[4] >= 0) mask += 2;
  if (conn[5] >= 0) mask += 4;

  switch (mask)
  {
    case  0: // 000
      tri[nt][0] = conn[0];
      tri[nt][1] = conn[1];
      tri[nt][2] = conn[2];
      nt++;
      break;
    case  1: // 001
      tri[nt][0] = conn[0];
      tri[nt][1] = conn[3];
      tri[nt][2] = conn[2];
      nt++;
      tri[nt][0] = conn[3];
      tri[nt][1] = conn[1];
      tri[nt][2] = conn[2];
      nt++;
      break;
    case  2: // 010
      tri[nt][0] = conn[0];
      tri[nt][1] = conn[1];
      tri[nt][2] = conn[4];
      nt++;
      tri[nt][0] = conn[0];
      tri[nt][1] = conn[4];
      tri[nt][2] = conn[2];
      nt++;
      break;
    case  3: // 011
      tri[nt][0] = conn[3];
      tri[nt][1] = conn[1];
      tri[nt][2] = conn[4];
      nt++;
      asp1 =          triangle_aspect_ratio(nodep[conn[0]],nodep[conn[3]],nodep[conn[4]]);
      asp1 = MAX(asp1,triangle_aspect_ratio(nodep[conn[0]],nodep[conn[4]],nodep[conn[2]]));
      asp2 =          triangle_aspect_ratio(nodep[conn[0]],nodep[conn[3]],nodep[conn[2]]);
      asp2 = MAX(asp2,triangle_aspect_ratio(nodep[conn[3]],nodep[conn[4]],nodep[conn[2]]));
      if (asp1 < asp2)
      {
        tri[nt][0] = conn[0];
        tri[nt][1] = conn[3];
        tri[nt][2] = conn[4];
        nt++;
        tri[nt][0] = conn[0];
        tri[nt][1] = conn[4];
        tri[nt][2] = conn[2];
        nt++;
      } else
      {
        tri[nt][0] = conn[0];
        tri[nt][1] = conn[3];
        tri[nt][2] = conn[2];
        nt++;
        tri[nt][0] = conn[3];
        tri[nt][1] = conn[4];
        tri[nt][2] = conn[2];
        nt++;
      }
      break;
    case  4: // 100
      tri[nt][0] = conn[0];
      tri[nt][1] = conn[1];
      tri[nt][2] = conn[5];
      nt++;
      tri[nt][0] = conn[5];
      tri[nt][1] = conn[1];
      tri[nt][2] = conn[2];
      nt++;
      break;
    case  5: // 101
      tri[nt][0] = conn[0];
      tri[nt][1] = conn[3];
      tri[nt][2] = conn[5];
      nt++;
      asp1 =          triangle_aspect_ratio(nodep[conn[3]],nodep[conn[1]],nodep[conn[2]]);
      asp1 = MAX(asp1,triangle_aspect_ratio(nodep[conn[3]],nodep[conn[2]],nodep[conn[5]]));
      asp2 =          triangle_aspect_ratio(nodep[conn[3]],nodep[conn[1]],nodep[conn[5]]);
      asp2 = MAX(asp2,triangle_aspect_ratio(nodep[conn[1]],nodep[conn[2]],nodep[conn[5]]));
      if (asp1 < asp2)
      {
        tri[nt][0] = conn[3];
        tri[nt][1] = conn[1];
        tri[nt][2] = conn[2];
        nt++;
        tri[nt][0] = conn[3];
        tri[nt][1] = conn[2];
        tri[nt][2] = conn[5];
        nt++;
      } else
      {
        tri[nt][0] = conn[3];
        tri[nt][1] = conn[1];
        tri[nt][2] = conn[5];
        nt++;
        tri[nt][0] = conn[1];
        tri[nt][1] = conn[2];
        tri[nt][2] = conn[5];
        nt++;
      }
      break;
    case  6: // 110
      tri[nt][0] = conn[5];
      tri[nt][1] = conn[4];
      tri[nt][2] = conn[2];
      nt++;
      asp1 =          triangle_aspect_ratio(nodep[conn[0]],nodep[conn[1]],nodep[conn[4]]);
      asp1 = MAX(asp1,triangle_aspect_ratio(nodep[conn[0]],nodep[conn[4]],nodep[conn[5]]));
      asp2 =          triangle_aspect_ratio(nodep[conn[0]],nodep[conn[1]],nodep[conn[5]]);
      asp2 = MAX(asp2,triangle_aspect_ratio(nodep[conn[1]],nodep[conn[4]],nodep[conn[5]]));
      if (asp1 < asp2)
      {
        tri[nt][0] = conn[0];
        tri[nt][1] = conn[1];
        tri[nt][2] = conn[4];
        nt++;
        tri[nt][0] = conn[0];
        tri[nt][1] = conn[4];
        tri[nt][2] = conn[5];
        nt++;
      } else
      {
        tri[nt][0] = conn[0];
        tri[nt][1] = conn[1];
        tri[nt][2] = conn[5];
        nt++;
        tri[nt][0] = conn[1];
        tri[nt][1] = conn[4];
        tri[nt][2] = conn[5];
        nt++;
      }
      break;
    case  7: // 111
      tri[nt][0] = conn[0];
      tri[nt][1] = conn[3];
      tri[nt][2] = conn[5];
      nt++;
      tri[nt][0] = conn[3];
      tri[nt][1] = conn[1];
      tri[nt][2] = conn[4];
      nt++;
      tri[nt][0] = conn[5];
      tri[nt][1] = conn[4];
      tri[nt][2] = conn[2];
      nt++;
      tri[nt][0] = conn[3];
      tri[nt][1] = conn[4];
      tri[nt][2] = conn[5];
      nt++;
      break;
    default:
      printf("\nREFINE_TRI: Tri marking not implemented, mask = %d",mask);
      fflush(stdout);
      abort();
      break;
  }

  return(nt);

}
