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
  ptiters = range(20,220,20)
  cfls = range(50,550,50)

  for ptiter in ptiters:
      for cfl in cfls:
	  mysystem('cp /home/vbetro/online_edu/images/order1/mach1.25/alpha-01/RMSEulerIter_1_-01_1.25_%03d_%03d_00.png /home/vbetro/online_edu/images/order2/mach1.25/alpha-01/RMSEulerIter_2_-01_1.25_%03d_%03d_00.png'%(ptiter,cfl,ptiter,cfl));
  for mach in machs:
      for alpha in alphas:
	  mysystem('cp /home/vbetro/online_edu/images/order1/mach%0.2f/alpha%+03d/CPEuler_3_1_%+03d_%0.2f_00.png /home/vbetro/online_edu/images/order2/mach%0.2f/alpha%+03d/CPEuler_3_2_%+03d_%0.2f_00.png'%(mach,alpha,alpha,mach,mach,alpha,alpha,mach));
	  mysystem('cp /home/vbetro/online_edu/images/order1/mach%0.2f/alpha%+03d/CPEuler_3_1_%+03d_%0.2f_01.png /home/vbetro/online_edu/images/order2/mach%0.2f/alpha%+03d/CPEuler_3_2_%+03d_%0.2f_01.png'%(mach,alpha,alpha,mach,mach,alpha,alpha,mach));	  
	  for ptiter in ptiters:
              for cfl in cfls:
		  if alpha == 0:
                      mysystem('cp /home/vbetro/online_edu/images/order1/mach%0.2f/alpha%+03d/RMSEulerIter_1_%+03d_%0.2f_%03d_%03d_00.png /home/vbetro/online_edu/images/order2/mach%0.2f/alpha%+03d/RMSEulerIter_2_%+03d_%0.2f_%03d_%03d_00.png'%(mach,alpha,alpha,mach,ptiter,cfl,mach,alpha,alpha,mach,ptiter,cfl));
                  else:
                      mysystem('cp /home/vbetro/online_edu/images/order2/mach%0.2f/alpha-%02d/RMSEulerIter_2_-%02d_%0.2f_%03d_%03d_00.png /home/vbetro/online_edu/images/order2/mach%0.2f/alpha%+03d/RMSEulerIter_2_%+03d_%0.2f_%03d_%03d_00.png'%(mach,alpha,alpha,mach,ptiter,cfl,mach,alpha,alpha,mach,ptiter,cfl));
		      
  
if __name__ == "__main__":
   main()
