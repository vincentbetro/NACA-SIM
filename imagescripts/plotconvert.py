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
  #now, we need to recursively move everybody back
  for order in orders:
      for mach in machs:
          for alpha in alphas:
                  mysystem('convert -auto-orient -density 100 -rotate 90 -resize 640x480 /ibrix-scr/vbetro/imagescripts/CPEuler_3_%d_%+03d_%0.2f_200_050_01.ps /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/CPEuler_3_%d_%+03d_%0.2f_02.png'%(order,alpha,mach,order,mach,alpha,order,alpha,mach));
                  mysystem('mv /ibrix-scr/vbetro/imagescripts/CPEuler_3_%d_%+03d_%0.2f_200_050_01.ps /ibrix-scr/vbetro/psimages/CPEuler_3_%d_%+03d_%0.2f_02.ps'%(order,alpha,mach,order,alpha,mach));

if __name__ == "__main__":
   main()
