#!/usr/bin/env python

import sys
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

# input xyz file
xyzfile = "output.traj.xyz"
infile = open(xyzfile, 'r')

# read input and store in the array xyz
xyz = []
_CONTINUE_ = True
while(_CONTINUE_):
    for i in range(3):
        line = infile.readline()
        if i==0:
            if(len(line.strip())>0):
                if(int(line)!=1):
                    sys.exit('plotTraj.py only works for single prmry aglmrts')
        if i==2:
            tmp = line.split()
            if len(tmp)>0:
                xyz.append([float(tmp[i]) for i in [1,2,3]])
        if len(line.strip())==0:
            _CONTINUE_ = False

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = Axes3D(fig) #fig.gca(projection='3d')

x = np.array(xyz)
ax.plot(x[:,0],x[:,1],x[:,2])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

#ax.set_xlim3d(0, 100)
#ax.set_ylim3d(0, 100)
#ax.set_zlim3d(0, 100)


plt.show()
