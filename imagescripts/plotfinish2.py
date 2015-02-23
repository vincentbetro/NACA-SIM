#!/usr/bin/env python

import os,sys
import glob

def mysystem(s):
    print(s)
    retval = os.system(s)
    return retval

def main():
  alphas = range(0,4)
  machs = [0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25]
  for mach in machs:
      for alpha in alphas:
	  mysystem('cp /home/vbetro/online_edu/images/order1/mach%0.2f/alpha%+03d/CPEuler_3_1_%+03d_%0.2f_00.png /home/vbetro/online_edu/images/order2/mach%0.2f/alpha%+03d/CPEuler_3_2_%+03d_%0.2f_00.png'%(mach,alpha,alpha,mach,mach,alpha,alpha,mach));
	  mysystem('cp /home/vbetro/online_edu/images/order1/mach%0.2f/alpha%+03d/CPEuler_3_1_%+03d_%0.2f_01.png /home/vbetro/online_edu/images/order2/mach%0.2f/alpha%+03d/CPEuler_3_2_%+03d_%0.2f_01.png'%(mach,alpha,alpha,mach,mach,alpha,alpha,mach));	  
		      
  
if __name__ == "__main__":
   main()
