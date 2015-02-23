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

  for order in orders:
      for alpha in alphas:
          for mach in machs:
                  mysystem("gracebat /ibrix-scr/vbetro/order%d/mach%0.2f/alpha%+03d/CPEuler_2_%d_%+03d_%0.2f_200_050_01.dat /ibrix-scr/vbetro/order%d/mach%0.2f/alpha%+03d/CPEuler_3_%d_%+03d_%0.2f_200_050_01.dat -p /ibrix-scr/vbetro/imagescripts/Cp.params"%(order,mach,alpha,order,alpha,mach,order,mach,alpha,order,alpha,mach));

if __name__ == "__main__":
    main()
