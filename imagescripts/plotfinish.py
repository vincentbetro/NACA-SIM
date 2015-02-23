#!/usr/bin/env python

import os,sys
import glob

def mysystem(s):
    print(s)
    retval = os.system(s)
    return retval

def main():
  alphas = range(-8,9)
  machs = [0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25]
  orders = [1,2]

  for mach in machs:
      for alpha in alphas:
          for order in orders:
              mysystem('mv /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/CPEuler_%d_%+03d_%0.2f_00.png /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/CPEuler_3_%d_%+03d_%0.2f_00.png'%(order,mach,alpha,order,alpha,mach,order,mach,alpha,order,alpha,mach));
              mysystem('mv /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/CPEuler_%d_%+03d_%0.2f_01.png /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/CPEuler_3_%d_%+03d_%0.2f_01.png'%(order,mach,alpha,order,alpha,mach,order,mach,alpha,order,alpha,mach));
		      
  
if __name__ == "__main__":
   main()
