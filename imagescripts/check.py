#!/usr/bin/env python

import os,sys
import glob

def mysystem(s):
    print(s)
    retval = os.system(s)
    return retval

def main():
    mysystem('ls -l | grep -c png');

if __name__ == "__main__":
   main()
