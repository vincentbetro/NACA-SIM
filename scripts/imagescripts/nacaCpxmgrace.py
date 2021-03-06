#!/usr/bin/env python

import os,sys
import glob

def mysystem(s):
    print(s)
    retval = os.system(s)
    return retval

def main():
  alphas = range(-10,11)
  orders = [1,2]
  machs = [0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25]
  #allows us to get whole range, excluding last number, and inc by third value
  cfls = range(50,550,50)
  ptiters = range(20,220,20)

  for order in orders:
      for alpha in alphas:
          for mach in machs:
              for cfl in cfls:
                  for ptiter in ptiters:
                      mysystem("gracebat CPEuler_2_%d_%+03d_%0.2f_%03d_%03d_00.dat CPEuler_2_%d_%+03d_%0.2f_%03d_%03d_00.dat -p Cp.params"%(order,alpha,mach,ptiter,cfl,order,alpha,mach,ptiter,cfl));

if __name__ == "__main__":
    main()
