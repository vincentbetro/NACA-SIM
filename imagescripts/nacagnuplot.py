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
  refinelevels = [2]

  for order in orders:
      for alpha in alphas:
          for mach in machs:
              for refinelevel in refinelevels:
                  mysystem("./Winslow.sh ../meshes/naca0012_%d_%0.2f_%+03d_%02d"%(order,mach,alpha,refinelevel));

if __name__ == "__main__":
    main()
