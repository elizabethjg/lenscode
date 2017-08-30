import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import os

""" Ejecuta sextractor. Calcula el seeing y el punto de saturacion. Edita estos parametros
en el archivo de configuracion y vuelve a correr sextractor para un threshold mas bajo. 

ARCHIVO DE ENTRADA: archivo_entrada

del cual lee, numero (parametro), nombre de la imagen, filtro, pixsize,
gain, mag_zeropoint, magmax, magmin, fwhmmax (los ultimos tres parametros son para el 
calculo del seeing)

Finalmente ejecuta sextractor y genera el archivo de input para la clasificacion de las fuentes
(2selectstars.py)

"""

archivo_entrada='image_parameters_filter_r.in'


def seeing(imagen,pixsize,corrida,filtro,gain,mag_zeropoint,magmax,magmin,fwhmmax,plot):

	

	
	
	#run sextrator

	sexfile='first'+filtro+str(corrida)+'.sex' #name of configuration file of SExtractor
	salida_sex='run1-'+filtro+str(corrida)+'.cat' #exit catalogue from SExtractor

	#print '----------------------------------------------------------'
	#print '           WRITING SExtractor CONFIGURATION FILE          '
	#print '=========================================================='
	
	os.system('rm '+sexfile)
	os.system('rm '+salida_sex)	
	f1=open(sexfile,'w')
	f1.write('#-------------------------------- Catalog ------------------------------------\n')
	f1.write('CATALOG_NAME     '+salida_sex+'   # name of the output catalog\n')
	f1.write('CATALOG_TYPE     ASCII_HEAD     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,\n')
	f1.write('                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC\n')
	f1.write('PARAMETERS_NAME  first.param     # name of the file containing catalog contents \n')
	f1.write('#------------------------------- Extraction ----------------------------------\n')
	f1.write('DETECT_TYPE      CCD            # CCD (linear) or PHOTO (with gamma correction)\n')
	f1.write('DETECT_MINAREA   5.0              # minimum number of pixels above threshold\n')
	f1.write('DETECT_THRESH    5.0            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n')
	f1.write('ANALYSIS_THRESH  2.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n')
	f1.write('FILTER           Y              # apply filter for detection (Y or N)?\n')
	f1.write('FILTER_NAME      gauss_3.0_5x5.conv # name of the file containing the filter\n')
	f1.write('DEBLEND_NTHRESH  32            # Number of deblending sub-thresholds\n')
	f1.write('DEBLEND_MINCONT  0.005          # Minimum contrast parameter for deblending\n')
	f1.write('CLEAN            Y              # Clean spurious detections? (Y or N)?\n')
	f1.write('CLEAN_PARAM      1.0            # Cleaning efficiency\n')
	f1.write('MASK_TYPE        CORRECT        # type of detection MASKing: can be one of\n')
	f1.write('                                # NONE, BLANK or CORRECT\n')
	f1.write('#------------------------------ Photometry -----------------------------------\n')
	f1.write('PHOT_APERTURES   5             # MAG_APER aperture diameter(s) in pixels\n')
	f1.write('PHOT_AUTOPARAMS  2.5, 3.5       # MAG_AUTO parameters: <Kron_fact>,<min_radius>\n')
	f1.write('PHOT_PETROPARAMS 2.0, 3.5       # MAG_PETRO parameters: <Petrosian_fact>,\n')
	f1.write('                                # <min_radius>\n')
	f1.write('SATUR_LEVEL      50000.0        # level (in ADUs) at which arises saturation\n')
	f1.write('MAG_ZEROPOINT    '+str(mag_zeropoint)+'    # magnitude zero-point\n')
	f1.write('MAG_GAMMA        4.0            # gamma of emulsion (for photographic scans)\n')
	f1.write('GAIN             '+str(gain)+'            # detector gain in e-/ADU\n')
	f1.write('PIXEL_SCALE      '+str(pixsize)+'         # size of pixel in arcsec (0=use FITS WCS info)\n')
	f1.write('#------------------------- Star/Galaxy Separation ----------------------------\n')
	f1.write('SEEING_FWHM      1.0            # stellar FWHM in arcsec\n')
	f1.write('STARNNW_NAME     default.nnw    # Neural-Network_Weight table filename\n')
	f1.write('#------------------------------ Background -----------------------------------\n')
	f1.write('BACK_SIZE        36             # Background mesh: <size> or <width>,<height>\n')
	f1.write('BACK_FILTERSIZE  8              # Background filter: <size> or <width>,<height>\n')
	f1.write('BACKPHOTO_TYPE   GLOBAL         # can be GLOBAL or LOCAL\n')
	f1.write('#------------------------------ Check Image ----------------------------------\n')
	f1.write('CHECKIMAGE_TYPE   APERTURES     # can be NONE, BACKGROUND, BACKGROUND_RMS,\n')
	f1.write('                                # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,\n')
	f1.write('                                # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,\n')
	f1.write('                                # or APERTURES\n')
	f1.write('CHECKIMAGE_NAME  07.fits         # Filename for the check-image\n')
	f1.write('#--------------------- Memory (change with caution!) -------------------------\n')
	f1.write('MEMORY_OBJSTACK  3000           # number of objects in stack\n')
	f1.write('MEMORY_PIXSTACK  300000         # number of pixels in stack\n')
	f1.write('MEMORY_BUFSIZE   1024           # number of lines in buffer\n')
	f1.write('#----------------------------- Miscellaneous ---------------------------------\n')
	f1.write('VERBOSE_TYPE     NORMAL         # can be QUIET, NORMAL or FULL\n')
	f1.write('WRITE_XML        N              # Write XML file (Y/N)?\n')
	f1.write('XML_NAME         sex.xml        # Filename for XML output\n')
	f1.close()




	callsex='sex '+imagen+' -c '+sexfile+' > sex_output'
		
	os.system(callsex)

	#print imagen
	#wait=raw_input('chequear y presionar ENTER')
	cat = np.loadtxt(salida_sex, comments='#') #reads the catalogue of the first run of sextractor
	
	#print '----------------------------------------------------------'
	#print '              COMPUTE SATURATION LEVEL                    '
	#print '=========================================================='


	#compute the saturation level
	FLUXMAX = cat[:,8]
	SATUR=0.8*FLUXMAX.max()
	#print 'SATUR_LEVEL =',SATUR
	
	#print '----------------------------------------------------------'
	#print '                COMPUTE THE SEEING                        '
	#print '=========================================================='


	#now, starts to compute the seeing
	
	FWHM = cat[:,5]
	FLAG = cat[:,11]
	MAGBEST = cat[:,6]
	nobj=len(MAGBEST)
	fwhm=np.zeros((nobj),float)
	mag=np.zeros((nobj),float)
		
	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		plt.plot(FWHM,MAGBEST, 'b.')
		plt.ylabel('MAG_BEST')
		plt.xlabel('FWHM')
		plt.show()
	
	
	j=0
	


	#make the cuts in magnitude to get mostly the stars
	for i in range(nobj):
		if MAGBEST[i] < magmax and MAGBEST[i] > magmin and FWHM[i] < fwhmmax:
			fwhm[j]=FWHM[i]
			mag[j]=MAGBEST[i]
			j=j+1
	
	fwhm=fwhm[:j]
	mag=mag[:j]	
	
	
	#and get the maximun value of the FWHM distribution
	
	
	step=0.05
	intervalos=int((fwhm.max()-fwhm.min())/step)
	interval1=fwhm.min()
	interval2=interval1+step
	maximo=0
	moda=90.
	for j in range(intervalos):
		contar=0
		for x in fwhm:
			if x > interval1 and x < interval2:
				contar=contar+1
			if contar > maximo:
				maximo=contar
				moda=(interval1+interval2)/2.
		interval1=interval1+step
		interval2=interval2+step
	
	
	#print 'SEEING in pix',moda
	#print 'SEEING in arcsec',moda*pixsize
	
	Seeing_file=np.round(moda)
	
	if Seeing_file < 3.:
		filter_file='gauss_2.0_5x5.conv'
	elif Seeing_file == 3.:
		filter_file='gauss_3.0_7x7.conv'
	elif Seeing_file == 4.:
		filter_file='gauss_4.0_7x7.conv'
	else:
		filter_file='gauss_5.0_9x9.conv'
		
	seeing=moda*pixsize
	print ' '
	print ' ------------ seeing: ',seeing
	print ' ------------ imagen: ',imagen
	print ' '
	puntox=np.array([moda],float)
	puntoy=np.array([18.],float)

	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		
		puntox=np.array([moda],float)
		puntoy=np.array([18.],float)
		plt.plot(FWHM,MAGBEST, 'b.',puntox,puntoy,'ro')
		plt.axvline(fwhmmax, color = 'r')
		plt.axhline(magmax, color = 'r')
		plt.axhline(magmin, color = 'r')
		plt.ylabel('MAG_BEST')
		plt.xlabel('FWHM')
		plt.show()
		
	
	
	#~ print '----------------------------------------------------------'
	#~ print '           WRITING SExtractor CONFIGURATION FILE          '
	#~ print '=========================================================='
	
	sexfile='second'+filtro+str(corrida)+'.sex' #name of configuration file of SExtractor
	salida_sex='run2-'+filtro+str(corrida)+'.cat' #exit catalogue from SExtractor

	
	os.system('rm '+salida_sex)
	os.system('rm '+sexfile)
	#~ print corrida
	
	f1=open(sexfile,'w')
	f1.write('#-------------------------------- Catalog ------------------------------------\n')
	f1.write('CATALOG_NAME     '+salida_sex+'   # name of the output catalog\n')
	f1.write('CATALOG_TYPE     ASCII_HEAD     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,\n')
	f1.write('                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC\n')
	f1.write('PARAMETERS_NAME  second.param     # name of the file containing catalog contents \n')
	f1.write('#------------------------------- Extraction ----------------------------------\n')
	f1.write('DETECT_TYPE      CCD            # CCD (linear) or PHOTO (with gamma correction)\n')
	f1.write('DETECT_MINAREA   5.0              # minimum number of pixels above threshold\n')
	f1.write('DETECT_THRESH    1.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n')
	f1.write('ANALYSIS_THRESH  2.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n')
	f1.write('FILTER           Y              # apply filter for detection (Y or N)?\n')
	f1.write('FILTER_NAME      '+filter_file+' # name of the file containing the filter\n')
	f1.write('DEBLEND_NTHRESH  60            # Number of deblending sub-thresholds\n')
	f1.write('DEBLEND_MINCONT  0.001          # Minimum contrast parameter for deblending\n')
	f1.write('CLEAN            Y              # Clean spurious detections? (Y or N)?\n')
	f1.write('CLEAN_PARAM      1.0            # Cleaning efficiency\n')
	f1.write('MASK_TYPE        CORRECT        # type of detection MASKing: can be one of\n')
	f1.write('                                # NONE, BLANK or CORRECT\n')
	f1.write('#------------------------------ Photometry -----------------------------------\n')
	f1.write('PHOT_APERTURES   5             # MAG_APER aperture diameter(s) in pixels\n')
	f1.write('PHOT_AUTOPARAMS  2.5, 3.5       # MAG_AUTO parameters: <Kron_fact>,<min_radius>\n')
	f1.write('PHOT_PETROPARAMS 2.0, 3.5       # MAG_PETRO parameters: <Petrosian_fact>,\n')
	f1.write('                                # <min_radius>\n')
	f1.write('SATUR_LEVEL      '+str(int(SATUR))+'        # level (in ADUs) at which arises saturation\n')
	f1.write('MAG_ZEROPOINT    '+str(mag_zeropoint)+'    # magnitude zero-point\n')
	f1.write('MAG_GAMMA        4.0            # gamma of emulsion (for photographic scans)\n')
	f1.write('GAIN             '+str(gain)+'            # detector gain in e-/ADU\n')
	f1.write('PIXEL_SCALE      '+str(pixsize)+'         # size of pixel in arcsec (0=use FITS WCS info)\n')
	f1.write('#------------------------- Star/Galaxy Separation ----------------------------\n')
	f1.write('SEEING_FWHM      '+str('%.3f' % seeing)+'            # stellar FWHM in arcsec\n')
	f1.write('STARNNW_NAME     default.nnw    # Neural-Network_Weight table filename\n')
	f1.write('#------------------------------ Background -----------------------------------\n')
	f1.write('BACK_SIZE        36             # Background mesh: <size> or <width>,<height>\n')
	f1.write('BACK_FILTERSIZE  8              # Background filter: <size> or <width>,<height>\n')
	f1.write('BACKPHOTO_TYPE   GLOBAL         # can be GLOBAL or LOCAL\n')
	f1.write('#------------------------------ Check Image ----------------------------------\n')
	f1.write('CHECKIMAGE_TYPE   APERTURES     # can be NONE, BACKGROUND, BACKGROUND_RMS,\n')
	f1.write('                                # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,\n')
	f1.write('                                # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,\n')
	f1.write('                                # or APERTURES\n')
	f1.write('CHECKIMAGE_NAME  07.fits         # Filename for the check-image\n')
	f1.write('#--------------------- Memory (change with caution!) -------------------------\n')
	f1.write('MEMORY_OBJSTACK  3000           # number of objects in stack\n')
	f1.write('MEMORY_PIXSTACK  300000         # number of pixels in stack\n')
	f1.write('MEMORY_BUFSIZE   1024           # number of lines in buffer\n')
	f1.write('#----------------------------- Miscellaneous ---------------------------------\n')
	f1.write('VERBOSE_TYPE     NORMAL         # can be QUIET, NORMAL or FULL\n')
	f1.write('WRITE_XML        N              # Write XML file (Y/N)?\n')
	f1.write('XML_NAME         sex.xml        # Filename for XML output\n')
	f1.close()
	
	corrida_sex='sex '+imagen+' -c '+sexfile
	os.system(corrida_sex)
	del(cat)
	del(FWHM)
	del(FLAG)
	del(MAGBEST)
	del(nobj)
	del(fwhm)
	del(mag)

	print '----------------------------------------------------------------'
	print '                      RUNNING SExtractor (second time)          '
	print '================================================================'

	
	
	os.system('rm input_2select'+filtro+'.in')
	f2=open('input_2select'+filtro+'.in','w')

	f2.write('#image name \n')
	f2.write(imagen+'\n')
	f2.write('#sextractor catalogue \n')
	f2.write(salida_sex+'\n')
	f2.write('#seeing en pixs \n')
	f2.write(str(moda))
	f2.close()
	
	
	
	
	return seeing



# READ INPUT FILE
read=np.loadtxt(archivo_entrada,comments='#',dtype='str')

corrida=read[0].astype(int)
imagen=read[1]
filtro=read[2]
pixsize=read[3].astype(float)
gain=read[4].astype(float)
mag_zeropoint=read[5].astype(float)
magmax=read[6].astype(float)
magmin=read[7].astype(float)
fwhmmax=read[8].astype(float)

print seeing (imagen,pixsize,corrida,filtro,gain,mag_zeropoint,magmax,magmin,fwhmmax,'si')
