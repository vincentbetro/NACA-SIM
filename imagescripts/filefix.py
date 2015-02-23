#!/usr/bin/env python

import os,sys
import glob

def mysystem(s):
    print(s)
    retval = os.system(s)
    return retval

def main():
  alphas = range(-10,-6)
  orders = [1,2]
  machs = [0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25]

  #collects all here so we can fix
  mysystem("mkdir order11");
  mysystem("mkdir order21");

  mysystem('find -name "RMSEulerIter_1_*.dat" -exec mv {} /tmp/vbetro/order11/ \;');
  mysystem('find -name "TecplotEuler_1_*.dat" -exec mv {} /tmp/vbetro/order11/ \;');
  mysystem('find -name "CPEuler_2_1_*.dat" -exec mv {} /tmp/vbetro/order11/ \;');
  mysystem('find -name "CPEuler_3_1_*.dat" -exec mv {} /tmp/vbetro/order11/ \;');
  mysystem('find -name "RMSEulerIter_2_*.dat" -exec mv {} /tmp/vbetro/order21/ \;');
  mysystem('find -name "TecplotEuler_2_*.dat" -exec mv {} /tmp/vbetro/order21/ \;');
  mysystem('find -name "CPEuler_2_2_*.dat" -exec mv {} /tmp/vbetro/order21/ \;');
  mysystem('find -name "CPEuler_3_2_*.dat" -exec mv {} /tmp/vbetro/order21/ \;');
 
  #now, we need to recursively move everybody back
  for order in orders:
      for mach in machs:
          for alpha in alphas:
              mysystem('find -name "RMSEulerIter_%d_%+03d_%0.2f_*.dat" -exec mv {} /tmp/vbetro/order%d/mach%0.2f/alpha%+03d \;'%(order,alpha,mach,order,mach,alpha));
              mysystem('find -name "TecplotEuler_%d_%+03d_%0.2f_*.dat" -exec mv {} /tmp/vbetro/order%d/mach%0.2f/alpha%+03d \;'%(order,alpha,mach,order,mach,alpha));
              mysystem('find -name "CPEuler_2_%d_%+03d_%0.2f_*.dat" -exec mv {} /tmp/vbetro/order%d/mach%0.2f/alpha%+03d \;'%(order,alpha,mach,order,mach,alpha));
              mysystem('find -name "CPEuler_3_%d_%+03d_%0.2f_*.dat" -exec mv {} /tmp/vbetro/order%d/mach%0.2f/alpha%+03d \;'%(order,alpha,mach,order,mach,alpha));

if __name__ == "__main__":
   main()
