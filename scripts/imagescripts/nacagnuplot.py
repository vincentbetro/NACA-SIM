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

  for order in orders:
      for alpha in alphas:
          for mach in machs:
              mysystem("cd /ibrix-scr/vbetro/order%d/mach%0.2f/alpha%+03d"%(order,mach,alpha))
              mysystem("./Winslow.sh naca0012_%d_%0.2f_%+03d_00"%(order,mach,alpha));

if __name__ == "__main__":
    main()
