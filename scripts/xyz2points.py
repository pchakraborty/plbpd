#!/usr/bin/env python
# xyz2points.py


"""
converts an xyz file to a data file whose each line is the particle coords
"""

# input xyz file
xyzfile = "traj.xyz"

# output
outputfile = "traj.dat"

# open files
infile = open(xyzfile, 'r')
outfile = open(outputfile, 'w')

# read input
_CONTINUE_ = True
while(_CONTINUE_):
    for i in range(3):
        line = infile.readline()
        if i==2:
            tmp = line.split()
            if len(tmp)>0:
                outfile.write(' %s %s %s\n' % (tmp[1], tmp[2], tmp[3]))
            else:
                _CONTINUE_ = False;
                          
infile.close()
outfile.close()
