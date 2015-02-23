#!/usr/bin/env python

#import sys libraries
import os,sys

#first, define system call
def mysystem(s):
    #want to be able to see statement executed
    print(s)
    os.system(s)

def main():
      mysystem("./gridmod -igrid sea_fighter_with_jets_half_inviscid_26.9M.crunch -ogrid sea_fighter_with_jets_half_inviscid_26.9M2.crunch -mirrorz -validate")
      mysystem("./gridglue sea_fighter_with_jets_half_inviscid_26.9M.crunch sea_fighter_with_jets_half_inviscid_26.9M2.crunch sea_fighter_with_jets_half_inviscid_26.9MWHOLE.crunch")
      mysystem("./Conv.LINUX64 Conv.inp")
      mysystem("./cp seafighter_whole.tri iso")
      mysystem("./cp seafighter_whole.tri anise")
      mysystem("./cd iso")
      mysystem("./FASTAR.LINUX64 FASTAR.inp > out.txt &")
      mysystem("./cd ../aniso")
      mysystem("./FASTAR.LINUX64 FASTAR.inp > out.txt &")
      mysystem("./cd /ibrix-scr/vbetro/ellipse/phugg/adaptH")
      mysystem("./Conv.LINUX64 Conv.inp")
      mysystem("./decomp -grid P_HUGG.crunch -case PHUGGELLH2 -np 8")
      mysystem("./ssh euler1d")
      mysystem("./ssh 0nJD-7884")
      mysystem("./pysub PHUGGELLH2 -qblades")

if __name__ == "__main__":
   main()
