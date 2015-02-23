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
              mysystem('cd /ibrix-scr/vbetro/order%d/mach%0.2f/alpha%+03d'%(order,mach,alpha));
              for ptiter in ptiters:
                  for cfl in cfls:
                      for refinelevel in refinelevels:
                          mysystem('grep "nan" *.dat > /ibrix-scr/vbetro/files_%d_%0.2f_%+03d.out'%(order,mach,alpha));

if __name__ == "__main__":
   main()
