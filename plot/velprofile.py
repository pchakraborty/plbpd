#!/usr/bin/env python

import sys
import argparse
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

def main():
    args = parseCommandLineArgs()
    print(args)
    
    a = h5.File(args.file)
    flowVelocity = a["FlowVelocity"]

    if (args.xsection == "y"):
        plotCrossSectionY(flowVelocity)
    elif (args.xsection == "x"):
        plotCrossSectionX(flowVelocity)
    else:
        raise Exception("Incorrect xsection [%s]" % xsection)
        
def plotCrossSectionY(flowVelocity):
    zdim, ydim, xdim, i = flowVelocity.shape
    stride = max(ydim//5, 1)
    for y in range(1, ydim-1, stride):
        # 1:zdim-1, 1:ydim-1, 1:xdim-1 excludes the buffer layer
        u = flowVelocity[1:zdim-1, y, 1:xdim-1, 0];
        v = flowVelocity[1:zdim-1, y, 1:xdim-1, 1];
        plt.quiver(u,v)
        plt.title("y = %d" % y)
        plt.show()

def plotCrossSectionX(flowVelocity):
    zdim, ydim, xdim, i = flowVelocity.shape
    stride = max(xdim//5, 1)
    for x in range(1, xdim-1, stride):
        # 1:zdim-1, 1:ydim-1, 1:xdim-1 excludes the buffer layer
        u = flowVelocity[1:zdim-1, 1:ydim-1, x, 0];
        v = flowVelocity[1:zdim-1, 1:ydim-1, x, 1];
        plt.quiver(u,v)
        plt.title("x = %d" % x)
        plt.show()
            
def parseCommandLineArgs():
    p = argparse.ArgumentParser(description="TODO")
    p.add_argument("file", help="HDF5 output file")
    p.add_argument("xsection", choices=["x", "y"], help="cross-section to plot")
    return p.parse_args()

if __name__=="__main__":
    main()
