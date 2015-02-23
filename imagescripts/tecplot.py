#!/usr/bin/env python

import os,sys
import glob

def mysystem(s):
    print(s)
    retval = os.system(s)
    return retval

def main():
  alphas = range(-8,9)
  orders = [2]
  machs = [0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25]
  refinelevels = [0,1,2]

  for order in orders:
      for alpha in alphas:
          for mach in machs:
              for refinelevel in refinelevels:
                  result = '/ibrix-scr/vbetro/order%d/mach%0.2f/alpha%+03d/TecplotEuler_%d_%+03d_%0.2f_200_500_%02d.dat'%(order,mach,alpha,order,alpha,mach,refinelevel);
                  if os.path.exists(result):
                      mysystem("ln -s %s temp.dat"%(result));
                      mysystem("tecplot -mesa /ibrix-scr/vbetro/imagescripts/naca0012.lay temp.dat -b -p /ibrix-scr/vbetro/imagescripts/export.mcr -y TecplotEuler_%d_%+03d_%0.2f_%02d.png"%(order,alpha,mach,refinelevel));
                      mysystem("rm -f temp.dat");

if __name__ == "__main__":
    main()
