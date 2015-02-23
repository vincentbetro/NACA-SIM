#!/usr/bin/env python

#import sys and global var libraries, as well as option parser to make command line args
import os,sys
import glob
from optparse import OptionParser

#first, define system call with default retval for debugging
def mysystem(s,defaultretval=0):
    #allows us to set effective debug flag
    global dryrun
    #want to be able to see statement executed
    print(s)
    #set for debugging
    retval = defaultretval
    if not dryrun:
        retval = os.system(s)
    return retval

#allows us to set up file names with all requisite info, including a default refinelevel of 0 and blank pre and post ext
#+03 gives signed three digit, auto filled integer
def format_name(case,order,mach,alpha,refinelevel=0,postext="",preext=""):
    return "%s_%d_%0.2f_%+03d_%02d%s%s"%(case,order,mach,alpha,refinelevel,preext,postext)

#writes file to rotate grids with multiline string
def write_smooth_file(fname,case,alpha):
    f = open(fname,"w")
    s = """%s_%+03d
0
1
0.25
100
1
1
%d
2 3
1.0
1.0e-14
"""%(case,alpha,alpha)
    f.write(s)

#writes opt smooth file
def write_smooth1_file(fname,case,order,mach,alpha,refinelevel):
    f = open(fname,"w")
    #allows us to concatenate carriage return
    f.write(format_name(case,order,mach,alpha,refinelevel)+"\n")
    s = """0
2
500
1
1
0
2 3
1.0e-6
"""
    f.write(s)

#writes subdiv file, always sets output one level higher
def write_subdiv_file(fname,case,order,mach,alpha,refinelevel):
    f = open(fname,"w")
    f.write(format_name(case,order,mach,alpha,refinelevel+1)+"\n")
    s = """1
3
5
0
1.0
1.0
1.0
5.0
5.0
5.0
2.0
2.0
2.0
1
"""
    f.write(s)

#writes euler file
def write_euler_file(fname,case,alpha,mach,order,cfl,ptiter,refinelevel,extension="",path=""):
    f = open(fname,"w")
    s = """%d
%0.2f
1.4
"""%(alpha,mach)
    f.write(s)
    f.write("%s%s"%(path,format_name(case,order,mach,alpha,refinelevel,".mesh",extension)) + "\n")
    s = """2 3
1
1.0
%d
100
%d
%d
10000
1.0
1.0e-15
"""%(cfl,order,ptiter)
    f.write(s)


def main():
  global dryrun
  parser = OptionParser()
  parser.add_option("--grids",action="store_true",dest="grids",default=False,help="Generates only initial grids at all alphas.  Parallel must be set to 0.")
  parser.add_option("--dryrun",action="store_true",dest="dryrun",default=False,help="Just print the commands; do not execute them.")
  parser.add_option("--case",dest="case",default="naca0012",help="Original meshfile name, without extension.")
  parser.add_option("--parallel",dest="parallel",default="0",help="Splits job into 21 separate jobs. Each must be given proc number from 1 to 21.  Zero may only be used for generating grids.")
  (options,args) = parser.parse_args()

  #sets global variable to allow retvals to reflect debug and not execute
  dryrun = options.dryrun

  #if we set parallel to 0, runs all on one
  #else, we need to split up parallel artifically (which could be done more automatically, but it is simple to do it this way too)
  if options.parallel == "0":
      alphas = range(-10,11)
  if options.parallel == "2":
      mach = 0.85
  orders = [1]
  alphas = [-6,-5,0,1,4,5,6,7,8]
  #allows us to get whole range, excluding last number, and inc by third value
  cfls = range(50,550,50)
  ptiters = range(20,220,20)

  #always do grid run separate
  if options.grids:
    for alpha in alphas:
      write_smooth_file("MYSMOOTH",options.case,alpha)
      mysystem("./SMOOTH %s.mesh %s.mesh < MYSMOOTH > stdo.out"%(options.case,options.case))

      for order in orders:
          for mach in machs:
              f1 = "%s_%+03d_01.dat"%(options.case,alpha)
              f2 = "/ibrix-scr/vbetro/meshes/%s"%format_name(options.case,order,mach,alpha,0,".dat")
              mysystem("cp %s %s"%(f1,f2))

              f1 = "%s_%+03d_01.mesh"%(options.case,alpha)
              f2 = "/ibrix-scr/vbetro/meshes/%s"%format_name(options.case,order,mach,alpha,0,".mesh")
              mysystem("cp %s %s"%(f1,f2))

    #now, remove all .dat and deprecated mesh files
    mysystem("rm -f *.dat *_01.mesh")
    sys.exit(1)

  #need to artifically set refinelevel
  refinelevel = 1

  #now, loop over all parameters and do all three adaptation runs for each
  for order in orders:
      for alpha in alphas:
          for cfl in cfls:
              for ptiter in ptiters:
                          write_euler_file("MYSOLVER%s"%options.parallel,options.case,alpha,mach,order,cfl,ptiter,refinelevel,"","/ibrix-scr/vbetro/meshes/")
                          result = mysystem("./EULER < MYSOLVER%s > stdo.out"%options.parallel)
                          #need to signify went fine without 1st 2nd switch
                          files = glob.glob("*_%d_%+03d_%0.2f_%03d_%03d_%02d.dat"%(order,alpha,mach,ptiter,cfl,refinelevel))
                          for f in files:
                              newf = f.replace(".dat","_00.dat")
                              mysystem("mv %s %s"%(f,newf))
                          #if we did not get results 2nd order, we do first then second and reappend name
                          if result==0 and order==2:
                              mysystem("rm -f *_%d_%+03d_%0.2f_%03d_%03d_%02d_00.dat"%(order,alpha,mach,ptiter,cfl,refinelevel))
                              write_euler_file("MYSOLVER%s"%options.parallel,options.case,alpha,mach,1,cfl,ptiter,refinelevel,"","/ibrix-scr/vbetro/meshes/")
                              mysystem("./EULER < MYSOLVER%s > stdo.out"%options.parallel)
                              mysystem("rm -f *_%d_%+03d_%0.2f_%03d_%03d_%02d.dat"%(order,alpha,mach,ptiter,cfl,refinelevel))
                              write_euler_file("MYSOLVER%s"%options.parallel,options.case,alpha,mach,order,cfl,ptiter,refinelevel,"_out")
                              result = mysystem("./EULER < MYSOLVER%s > stdo.out"%options.parallel)
                              files = glob.glob("*_%d_%+03d_%0.2f_%03d_%03d_%02d.dat"%(order,alpha,mach,ptiter,cfl,refinelevel))
                              for f in files:
                                  newf = f.replace(".dat","_12.dat")
                                  mysystem("mv %s %s"%(f,newf))

                          if result==0:
                              files = glob.glob("*_%d_%+03d_%0.2f_%03d_%03d_%02d*.dat"%(order,alpha,mach,ptiter,cfl,refinelevel))
                              for f in files:
                                  newf = f.replace(".dat","_nan.dat")
                                  mysystem("mv %s %s"%(f,newf))

                          if result==-1:
                              files = glob.glob("*_%d_%+03d_%0.2f_%03d_%03d_%02d*.dat"%(order,alpha,mach,ptiter,cfl,refinelevel))
                              for f in files:
                                  newf = f.replace(".dat","_uncvg.dat")
                                  mysystem("mv %s %s"%(f,newf))

                          d = "/tmp/vbetro/order%d/mach%0.2f/alpha%+03d"%(order,mach,alpha)
                          mysystem("mkdir -p " + d)
                          mysystem("mv *_%d_%+03d_%0.2f_%03d_%03d_%02d*.dat"%(order,alpha,mach,ptiter,cfl,refinelevel)  + d)

                          if result==1 and refinelevel < 2:
                              write_subdiv_file("MYSUBDIV%s"%options.parallel,options.case,order,mach,alpha,refinelevel)
                              fname = format_name(options.case,order,mach,alpha,refinelevel,".mesh","_out")
                              mysystem("./SMOOTH /ibrix-scr/vbetro/meshes/%s /ibrix-scr/vbetro/meshes/%s < MYSUBDIV%s > stdo.out"%(fname,fname,options.parallel))

                              write_smooth1_file("MYSMOOTH1%s"%options.parallel,options.case,order,mach,alpha,refinelevel+1)
                              fname = format_name(options.case,order,mach,alpha,refinelevel+1,".mesh") 
                              mysystem("./SMOOTH /ibrix-scr/vbetro/meshes/%s /ibrix-scr/vbetro/meshes/%s < MYSMOOTH1%s > stdo.out"%(fname,fname,options.parallel))

                              base =  format_name(options.case,order,mach,alpha,refinelevel+1) 
                              mysystem("mv %s_01.dat /ibrix-scr/vbetro/meshes/%s.dat"%(base,base))
                              mysystem("mv %s_01.mesh /ibrix-scr/vbetro/meshes/%s.mesh"%(base,base))

if __name__ == "__main__":
   main()
