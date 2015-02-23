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

  #now, we need to recursively move everybody back
  for order in orders:
      for mach in machs:
          for alpha in alphas:
              result = '/home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/CPEuler_3_%d_%+03d_%0.2f_00.png'%(order,mach,alpha,order,alpha,mach);
              result1 = '/home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/CPEuler_3_%d_%+03d_%0.2f_01.png'%(order,mach,alpha,order,alpha,mach);
              if os.path.exists(result) and not os.path.exists(result1):
                  mysystem('cp /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/CPEuler_3_%d_%+03d_%0.2f_00.png /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/CPEuler_%d_%+03d_%0.2f_01.png'%(order,mach,alpha,order,alpha,mach,order,mach,alpha,order,alpha,mach));
              if not os.path.exists(result) and os.path.exists(result1):
                  mysystem('cp /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/CPEuler_3_%d_%+03d_%0.2f_01.png /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/CPEuler_%d_%+03d_%0.2f_00.png'%(order,mach,alpha,order,alpha,mach,order,mach,alpha,order,alpha,mach));
              for ptiter in ptiters:
                  for cfl in cfls:
                      result = '/home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/RMSEulerIter_%d_%+03d_%0.2f_%03d_%03d_00.png'%(order,mach,alpha,order,alpha,mach,ptiter,cfl);
                      if os.path.exists(result):
                          mysystem('mv /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/RMSEulerIter_%d_%+03d_%0.2f_%03d_%03d_01.png /ibrix-scr/vbetro/extraRMS'%(order,mach,alpha,order,alpha,mach,ptiter,cfl));
                      if not os.path.exists(result):
                          mysystem('mv /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/RMSEulerIter_%d_%+03d_%0.2f_%03d_%03d_01.png /home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/RMSEulerIter_%d_%+03d_%0.2f_%03d_%03d_00.png'%(order,mach,alpha,order,alpha,mach,ptiter,cfl,order,mach,alpha,order,alpha,mach,ptiter,cfl));
  for order in orders:
      for mach in machs:
          for alpha in alphas:
              files = glob.glob('/home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d/*.png'%(order,mach,alpha));
              n = len(files);
              print('/home/vbetro/online_edu/images/order%d/mach%0.2f/alpha%+03d'%(order,mach,alpha));
              print('n = %d'%(n));

if __name__ == "__main__":
   main()
