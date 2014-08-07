#!/usr/bin/env python

import os
import sys

import pylab as pl

for num in range(0,21000,1000):
    paddednum = '%06d' % num
    cmd = '/home/pc/Research/codes/plbpd/scripts/plotMoments.py -f ' + \
        'output.moments.' + str(num) + ' -o vel' + paddednum + '.png' + \
        ' -y 0 -s 0.03'
    print cmd
    os.system(cmd)
    
convert2gif = 'convert -delay 10 -loop 1 vel*.png velprofile.gif'
print convert2gif
os.system(convert2gif)

