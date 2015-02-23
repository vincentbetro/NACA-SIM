#!/usr/bin/env python

import os,sys
import glob

def mysystem(s):
    print(s)
    retval = os.system(s)
    return retval

def main():
  alphas = range(-3,4)
  machs = [0.55,0.65,0.75,0.85,0.95]
  orders = [2]

  for mach in machs:
      for alpha in alphas:
          for order in orders:
              result = '/home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/TecplotEuler_%d_%+03d_%0.2f_02.png'%(order,mach,alpha,order,alpha,mach);
              if not os.path.exists(result):
	          mysystem('cp /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/TecplotEuler_%d_%+03d_%0.2f_01.png /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/TecplotEuler_%d_%+03d_%0.2f_02.png'%(order,mach,alpha,order,alpha,mach,order,mach,alpha,order,alpha,mach));
		      
  
if __name__ == "__main__":
   main()
