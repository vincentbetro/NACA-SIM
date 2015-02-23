#!/usr/bin/env python

#import sys libraries
import os,sys

#first, define system call
def mysystem(s):
    #want to be able to see statement executed
    print(s)
    os.system(s)

def main():
      mysystem("./adaptationruns1.py --case naca0012 --parallel 20 &")
      mysystem("./adaptationruns2.py --case naca0012 --parallel 16 &")
      mysystem("./adaptationruns5.py --case naca0012 --parallel 1 &")
      mysystem("./adaptationruns6.py --case naca0012 --parallel 2 &")

if __name__ == "__main__":
   main()
