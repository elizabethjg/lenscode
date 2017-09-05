import pyfits
import math
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

""" 

Realiza la clasificacion en forma interactiva de las fuentes de acuerdo
al catalogo de salida del sextractor.

Entrada: archivo_entrada

"""

archivo_entrada='input_2selectr.in'

def star_gx(salida_sex,plot,fwhm,mumin,mumax):
	
	
	print '######## CLASIFICANDO LOS OBJETOs #######'
	
	cat = np.loadtxt(salida_sex, comments='#') #reads the catalogue of the first run of sextractor

	#extract some nobjameters from the catalogue
	FWHM = cat[:,20]
	MAGBEST = cat[:,11]
	MUMAX = cat[:,23]
	CLASS=cat[:,24]
	FLAG=cat[:,25]
	nobj=len(MAGBEST)
	ALFA=cat[:,3]
	DELTA=cat[:,4]
	control2=np.zeros(4,float)
	comnobjar=0.0
	
	# STARS CANDIDATES
	mask_candidatas = (CLASS > 0.85)*(FWHM > fwhm-0.5)*(FWHM < fwhm+0.8)
	candidatas=cat[mask_candidatas]
	i=len(candidatas)	
	print 'n candidatas',i
	#~if plot in ('s', 'S', 'si', 'Si', 'SI'):
		#~print 'candidatas ',i#'m,n',m, n
		#~print 'seeing', fwhm
		#~plt.plot(MAGBEST,MUMAX, 'k.')
		#~plt.xlabel('MAG_BEST')
		#~plt.ylabel('MU_MAX')
		#~plt.show()
		#~
		#~plt.plot(MAGBEST,CLASS, 'k.')
		#~plt.xlabel('MAG_BEST')
		#~plt.ylabel('MU_MAX')
		#~plt.show()
	#~
		#~plt.plot(FWHM,MAGBEST, 'k.')
		#~plt.ylabel('MAGBEST')
		#~plt.xlabel('FWHM')
		#~plt.show()
    
	X=candidatas[:,11]
	Y=candidatas[:,23]
	stars_fit = lambda x,m,n: x*m+n
	popt, pcov = curve_fit(stars_fit, X, Y)
	m,n = popt	
	varx = np.array([X.min(), X.max()])
	vary = m*varx+n
	# x y of the line to do the selection
	#x=range(int(MAGBEST.max()-MAGBEST.min()))+MAGBEST.min()
	zero=MAGBEST.min()
	x=range(15)+zero
	y=m*x+n
	print 'm,n',m,n	
	ancho=0.4	#width in magnitudes 
		
	#print 'ancho',ancho
    #PARA SDSS
	#~ mumin=7.8 #float(raw_input('Ingrese mu minimo '))
	#~ mumax=20.5
    
	mu=m*MAGBEST+n
	mask_good = (MUMAX > mumin) * (FWHM > fwhm-0.5) * (FLAG < 4.0) # Quitamos falsas detecciones
	mask_stars = (MUMAX < mu+ancho) * (MUMAX < mumax) * (FWHM < fwhm+1.0) * (MUMAX > mu-ancho) * mask_good	# Seleccionamos las estrellas
	mask_gx = (~mask_stars) * (CLASS < 0.8) * ((MUMAX > mumax)+(MUMAX > mu-ancho)) *mask_good		# Seleccionamos las galaxias
	stars = cat[mask_stars]
	gx = cat[mask_gx]
	fake= cat[~mask_gx*~mask_stars]


	print 'cantidad de estrellas seleccionadas: ', len(stars)
	print 'cantidad de galaxias seleccionadas: ', len(gx)	
	
	starsmag=stars[:,11]
	starsmu=stars[:,23]
	starsfw=stars[:,20]
	
	gxmag=gx[:,11]
	gxmu=gx[:,23]
	gxfw=gx[:,20]
	
	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		print 'm,n',m, n
		#~plt.plot(MAGBEST,MUMAX, 'k.',starsmag,starsmu,'b.',gxmag, gxmu, 'g.',x,y,'r',varx,vary,'b')
		#~plt.xlabel('MAG_BEST')
		#~plt.ylabel('MU_MAX')
		#~plt.show()
	
		#~plt.plot(FWHM,MAGBEST, 'k.',starsfw,starsmag,'b*',gxfw,gxmag,'g.')
		#~plt.ylabel('MAGBEST')
		#~plt.xlabel('FWHM')
		#~plt.show()
#~
		y = np.linspace(mumin,mumax, 100)
		x=(y-n)/m
		x1 = ((y-n)/m)-ancho
		x2 = ((y-n)/m)+ancho
	
		mag3=((mumin-n)/m)+ancho
		mag4=((mumax-n)/m)-ancho
		x3 = np.linspace(MAGBEST.min(),mag3, 100)
		y3=np.empty(100)
		y3[:]=mumin

		x4 = np.linspace(mag4,MAGBEST.max(), 100)
		y4=np.empty(100)
		y4[:]=mumax

		rcParams['font.family'] = 'serif'
		rcParams['figure.figsize'] = 9., 8.0
		fig, ax = plt.subplots()
		majorFormatter = FormatStrFormatter('%.1f')

		plt.plot(fake[:,11],fake[:,23], 'kx')
		plt.plot(starsmag,starsmu,'b.',gxmag, gxmu, 'g.',x4,y4,'r',x3,y3,'r',x1,y,'r',x2,y,'r',alpha=0.6)
		plt.xlabel(u'$r^{\prime}$',fontsize=17)
		plt.ylabel('$\mu_{max}$ [mag arcsec$^{-2}$]',fontsize=17)
		#~ plt.axis([10,25,16,24])
		plt.subplots_adjust(left=0.1, right=0.85, top=0.85, bottom=0.1)
		plt.savefig('mu_mag_a2029.eps', format='eps',bbox_inches='tight')
		plt.show()	
		
		
		plt.plot(fake[:,20],fake[:,11], 'kx')
		plt.plot(starsfw,starsmag,'b.',gxfw,gxmag,'g.',alpha=0.6)
		#~ plt.axis([0,24,16,24])
		plt.ylabel('$r^{\prime}$',fontsize=17)
		plt.xlabel('FWHM (PIX)',fontsize=17)
		plt.subplots_adjust(left=0.1, right=0.85, top=0.85, bottom=0.1)
		plt.savefig('mag_fwhm_a2029.eps', format='eps',bbox_inches='tight')
		plt.show()
    #~
	return stars, gx



#READ CONFIGURATION FILE
read=np.loadtxt(archivo_entrada,comments='#',dtype='str')


image=read[0] #name of the image
salida_sex=read[1] #exit catalogue from SExtractor
fwhm=read[2].astype(float) #en pix




cat = np.loadtxt(salida_sex, comments='#') #reads the catalogue of the first run of sextractor

#extract some parameters from the catalogue
FWHM = cat[:,20]
MAGBEST = cat[:,11]
MUMAX = cat[:,23]
CLASS=cat[:,24]
FLAG=cat[:,25]
nobj=len(MAGBEST)
#some parameters



print '----------------------------------------------------------------'
print '                    START THE SELECTION                         '
print '================================================================'




plt.plot(MAGBEST,MUMAX, 'k.')
plt.xlabel('MAG_BEST')
plt.ylabel('MU_MAX')

plt.show()


mumin=float(raw_input('Ingrese mu minimo '))
mumax=float(raw_input('Ingrese mu maximo '))

stars, gx=star_gx(salida_sex,'si',fwhm,mumin,mumax)		

f1=open('stars.cat','w')

f1.write('26 ')
f1.write(str(len(stars)))
f1.write('\n')
np.savetxt(f1, stars, fmt='%6i   %8.3f   %7.3f   %10.7f	%10.7f	%7.3f	%7.3f	%5.1f	 %7.4f	 %7.4f  %5.1f   %8.4f	 %8.4f	 %8.4f	%8.4f	%8.4f	%8.4f	%8.4f	%8.4f 	%8i	%8.2f	%8.3f	%10.4f	  %8.4f	  %5.2f %2i')

f2=open('gx.cat','w')

f2.write('26 ')
f2.write(str(len(gx)))
f2.write('\n')

np.savetxt(f2, gx, fmt='%6i   %8.3f   %7.3f   %10.7f	%10.7f	%7.3f	%7.3f	%5.1f	 %7.4f	 %7.4f  %5.1f   %8.4f	 %8.4f	 %8.4f	%8.4f	%8.4f	%8.4f	%8.4f	%8.4f 	%8i	%8.2f	%8.3f	%10.4f	  %8.4f	  %5.2f %2i')


print 'END OF THE PROGRAM :)'

