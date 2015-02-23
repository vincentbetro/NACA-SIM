#!/usr/bin/env python

import os,sys
import glob

def mysystem(s):
    print(s)
    retval = os.system(s)
    return retval

def main():
  alphas = range(-8,9)
  orders = [1]
  machs = [0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25]
  ptiters = range(20,220,20)
  cfls = range(50,550,50)

  #now, we need to recursively move everybody back
  for order in orders:
      for mach in machs:
          for alpha in alphas:
              for ptiter in ptiters:
                  for cfl in cfls:
                      mysystem('find -name "RMSEulerIter_%d_%+03d_%0.2f_%03d_%03d.dat" -exec mv {} /ibrix-scr/vbetro/order%d/mach%0.2f/alpha%+03d/RMSEulerIter_%d_%+03d_%0.2f_%03d_%03d_01.dat \;'%(order,alpha,mach,ptiter,cfl,order,mach,alpha,order,alpha,mach,ptiter,cfl));
                      mysystem('find -name "TecplotEuler_%d_%+03d_%0.2f_%03d_%03d.dat" -exec mv {} /ibrix-scr/vbetro/order%d/mach%0.2f/alpha%+03d/TecplotEuler_%d_%+03d_%0.2f_%03d_%03d_01.dat \;'%(order,alpha,mach,ptiter,cfl,order,mach,alpha,order,alpha,mach,ptiter,cfl));
                      mysystem('find -name "CPEuler_2_%d_%+03d_%0.2f_%03d_%03d.dat" -exec mv {} /ibrix-scr/vbetro/order%d/mach%0.2f/alpha%+03d/CPEuler_2_%d_%+03d_%0.2f_%03d_%03d_01.dat \;'%(order,alpha,mach,ptiter,cfl,order,mach,alpha,order,alpha,mach,ptiter,cfl));
                      mysystem('find -name "CPEuler_3_%d_%+03d_%0.2f_%03d_%03d.dat" -exec mv {} /ibrix-scr/vbetro/order%d/mach%0.2f/alpha%+03d/CPEuler_3_%d_%+03d_%0.2f_%03d_%03d_01.dat \;'%(order,alpha,mach,ptiter,cfl,order,mach,alpha,order,alpha,mach,ptiter,cfl));

if __name__ == "__main__":
   main()
