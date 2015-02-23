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
  for order in orders:
      for mach in machs:
          for alpha in alphas:
              files = glob.glob('/Users/vincentbetro/Documents/school-sim/online_edu/images/order%d/mach%0.2f/alpha%+03d/*.png'%(order,mach,alpha));
              n = len(files);
              if n != 108:
                  print('/home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d'%(order,mach,alpha));
                  print('n = %d'%(n));

if __name__ == "__main__":
   main()
