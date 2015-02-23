#!/usr/bin/env python

import os,sys
import glob

def mysystem(s):
    print(s)
    retval = os.system(s)
    return retval

def main():
  alphas = range(-3,4)
  orders = [2]
  machs = [0.55,0.65,0.75,0.85,0.95]

  #now, we need to recursively move everybody back
  for order in orders:
      for mach in machs:
          for alpha in alphas:
              result = '/home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/naca0012_%d_%0.2f_%+03d_02.png'%(order,mach,alpha,order,mach,alpha);
              if not os.path.exists(result):   
                  mysystem('cp /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/naca0012_%d_%0.2f_%+03d_01.png /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/naca0012_%d_%0.2f_%+03d_02.png'%(order,mach,alpha,order,mach,alpha,order,mach,alpha,order,mach,alpha));

if __name__ == "__main__":
   main()
