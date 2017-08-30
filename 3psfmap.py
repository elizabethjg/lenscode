import math
import numpy as np
from pylab import *
import os
import pyfits
import sys
from selectstars import *
from psffile import *

f = open('stars.out', 'r')
ncol = int(f.read(2))
nstars = int(f.read(4))
f.close()
print ncol, nstars
allstars = np.loadtxt('stars.out', comments='#',skiprows=1) #reads the im2shape catalogue for the stars

f = open('gx.cat', 'r')
ncolcat = int(f.read(2))
ngxcat = int(f.read(7))
psfgx=np.zeros((ngxcat,8),float)
f.close()

starscat = np.loadtxt('stars.cat', comments='#',skiprows=1) #reads the stars from sextractor cat

print 'n gx',ngxcat
gxcat = np.loadtxt('gx.cat', comments='#',skiprows=1) #reads the stars from sextractor cat


goodstars=selectstars(allstars,'si')

out1=psffile(goodstars,starscat,5,'psfstars.dat')
out2=psffile(goodstars,gxcat,5,'psfgx.dat')

