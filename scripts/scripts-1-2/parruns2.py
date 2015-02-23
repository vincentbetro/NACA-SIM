#!/usr/bin/env python

#import sys libraries
import os,sys

#first, define system call
def mysystem(s):
    #want to be able to see statement executed
    print(s)
    os.system(s)

def main():
      mysystem("./adaptationruns2.py --case naca0012 --parallel 16 &")

if __name__ == "__main__":
   main()
