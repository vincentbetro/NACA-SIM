#!/usr/bin/env python

import os,sys
import glob

def mysystem(s):
    print(s)
    retval = os.system(s)
    return retval

def main():
  alphas = range(-8,9)
  orders = [1,2]
  machs = [0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25]
  ptiters = range(20,220,20)
  cfls = range(50,550,50)
  refinelevels = [0,1,2]

  #now, we need to recursively move everybody back
  for order in orders:
      for mach in machs:
          for alpha in alphas:
              for ptiter in ptiters:
                  for cfl in cfls:
                      for refinelevel in refinelevels:
                          #mysystem('mv CPEuler_3_%d_%+03d_%0.2f_%03d_%03d_%02d.ps /ibrix-scr/vbetro/images/order%d/mach%0.2f/alpha%+03d/CPEuler_3_%d_%+03d_%0.2f_%02d.ps'%(order,alpha,mach,ptiter,cfl,refinelevel,order,mach,alpha,order,alpha,mach,refinelevel));
                          #mysystem('mv TecplotEuler_%d_%+03d_%0.2f_%03d_%03d_%02d.png /ibrix-scr/vbetro/images/order%d/mach%0.2f/alpha%+03d/TecplotEuler_%d_%+03d_%0.2f_%02d.png'%(order,alpha,mach,ptiter,cfl,refinelevel,order,mach,alpha,order,alpha,mach,refinelevel));
                          mysystem('mv RMSEulerIter_%d_%+03d_%0.2f_%03d_%03d_%02d.ps /ibrix-scr/vbetro/images/order%d/mach%0.2f/alpha%+03d'%(order,alpha,mach,ptiter,cfl,refinelevel,order,mach,alpha));
                          #if refinelevel == 0:
                              #mysystem('mv /ibrix-scr/vbetro/meshes/naca0012_%d_%0.2f_%+03d_%02d.png /ibrix-scr/vbetro/images/meshes/naca0012_%+03d.png'%(order,mach,alpha,refinelevel,alpha));
                          #if refinelevel == 1:
                              #mysystem('mv /ibrix-scr/vbetro/meshes/naca0012_%d_%0.2f_%+03d_%02d.png /ibrix-scr/vbetro/images/order%d/mach%0.2f/alpha%+03d'%(order,mach,alpha,refinelevel,order,mach,alpha));
                          #if refinelevel == 2:
                              #mysystem('mv /ibrix-scr/vbetro/meshes/naca0012_%d_%0.2f_%+03d_%02d.png /ibrix-scr/vbetro/images/order%d/mach%0.2f/alpha%+03d'%(order,mach,alpha,refinelevel,order,mach,alpha));

if __name__ == "__main__":
   main()
