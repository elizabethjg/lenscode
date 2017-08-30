import pyfits
import math
import numpy as np
from pylab import *
import matplotlib.pyplot as plt

f = open('stars_correc.out', 'r')
ncol = int(f.read(2))
nstars = int(f.read(4))
f.close()
print ncol, nstars
allstars = np.loadtxt('stars_correc.out', comments='#',skiprows=1) #reads the im2shape catalogue for the stars
print allstars.shape

e1=allstars[:,18]
ab=allstars[:,10]
e2=allstars[:,19]
e=(e1**2+e2**2)**0.5
theta=(np.arctan2(e2,e1)/2.0)+(np.pi/2.0)
a=(ab*((1.0+e)/(1.0-e)))**0.5
posx= allstars[:,2]+allstars[:,6]
posy= allstars[:,3]+allstars[:,7]
ab1=allstars[:,10]


a1=a*cos(theta)
a2=a*sin(theta)

#make the check plots
plt.plot(e,ab,'k.')
plt.xlabel('e')
plt.ylabel('ab')
plt.axis([0,1,0,3])
plt.show()

plt.plot(e1,e2,'k.')
plt.xlabel('e1')
plt.ylabel('e2')
plt.title('Ellipticities after PSF correction')
plt.show()

posygraf=insert(posy,0,2000)
posxgraf=insert(posx,0,10.0)
a1graf=insert(a1,0,3.0)
a2graf=insert(a2,0,0.0)

quiver(posxgraf,posygraf,a1graf,a2graf,headwidth=1,headlength=0,scale_units='xy')
plt.title('Mayor axis map after PSF correction')
show()

quiver(posx,posy,e1,e2,headwidth=1,headlength=0,scale_units='xy')
plt.title('Ellipticity map after PSF correction')
show()

plt.hist(e1, 7, facecolor='b')
plt.title('e1 distribution after PSF correction')
plt.xlabel('e1')
show()

plt.hist(e2, 7, facecolor='b')
plt.title('e2 distribution after PSF correction')
plt.xlabel('e2')
show()

thetagrad=(theta*180.0)/np.pi
plt.hist(thetagrad, 7, facecolor='g')
plt.title('Theta distribution after PSF correction')
plt.xlabel('theta')
show()

plt.hist(a, 7, facecolor='g')
plt.title('a*b distribution after PSF correction')
plt.xlabel('a*b')
show()





print 'END OF THE PROGRAM :)'
