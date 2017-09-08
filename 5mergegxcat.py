import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import cosmolopy.distance as cd
import os
cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.7}

#EDIT THIS PART
pixsize=0.1454
z_cluster=0.501
cluster= 001
D_ang=cd.angular_diameter_distance(z_cluster, z0=0, **cosmo)
kpcscale=D_ang*(((1.0/3600.0)*np.pi)/180.0)*1000.0




#---------------------------------Reads cataloges----------------------------
f = open('gx.cat', 'r')
linecat=f.readline()
ncolcat = int(linecat[:2])
ngxcat = int(linecat[3:-1])
print ncolcat,ngxcat
gxcat = np.loadtxt(f, comments='#') #sextractor catalogue for galaxies
f.close()

cat2 = np.loadtxt('run2-i.cat', comments='#') #sextractor catalogue for all the objects in r band



f = open('gx.out', 'r') #im2shape catalog
linecat=f.readline()
ncol = int(linecat[:2])
ngxcat = int(linecat[3:-1])
gxim2 = np.loadtxt(f, comments='#') #im2shape out for galaxies
f.close()


#-----------------------EXTRACT AND COMPUTE SOME PARAMETERS--------------------------------

#parameters from sextractor's catalogues
idsex=gxcat[:,0]
x=gxcat[:,1]
y=gxcat[:,2]
alfa=gxcat[:,3]
delta=gxcat[:,4]
MAG1=gxcat[:,11]
err_MAG1=gxcat[:,12]
fwhm=gxcat[:,20]
mag2=cat2[:,11]
MUMAX=gxcat[:,23]
alfa2=cat2[:,3]
delta2=cat2[:,4]


#parameters from im2shape catalogue
idsex2=gxim2[:,1]
idim2=gxim2[:,0]
posx= gxim2[:,2]+gxim2[:,6]
posy= gxim2[:,3]+gxim2[:,7]
e1=gxim2[:,18]
e2=gxim2[:,19]
erre1=gxim2[:,34]
erre2=gxim2[:,35]
ab=gxim2[:,10]



ro_gx=ngxcat/(x.max()*y.max()*(pixsize**2)*((1.0/60.0)**2))
print 'densidad de gx = ', ro_gx, ' gx/arcmin**2'

e=(e1**2+e2**2)**0.5
theta=(np.arctan2(e2,e1)/2.0)+(np.pi/2.0)

print 'TOTAL DE GALAXIAS ',ngxcat
#Now is going to match sextractor catalogues to exctrat the r-magnitudes for the galaxies
MAG2=np.zeros(ngxcat,float)

difpos=np.zeros(ngxcat,float)
n=0

for i in range(ngxcat):
	difpos=((alfa2-alfa[i])**2+(delta2-delta[i])**2)**0.5
	MAG2[i]=mag2[argsort(difpos)[0]]	
	n=n+1
MAG2=MAG2[:n]



MAGBEST=MAG1



COLOR=MAG1-MAG2




print '------------------------------------------------------------------------'
print '                    COMPUTING THE MAX OF MAGBESTr                       '
print '========================================================================'


#compute the maximun of the r-magnitud distribution
#~ moda=stats.mode(MAGBEST)

#compute the maximun of the r-magnitud distribution
step=0.05

intervalos=int((MAGBEST.max()-MAGBEST.min())/step)
interval1=MAGBEST.min()
maximo=0

for j in range(intervalos):
	mask_mag=(MAGBEST>interval1)*(MAGBEST<(interval1+step))
	contar=len(MAGBEST[mask_mag])
	if contar > maximo:
		maximo=contar
		moda=(interval1+(interval1+step))/2.
	interval1=interval1+step



print 'maximo de la distribucion de magnitudes ',moda
plt.hist(MAGBEST,50,alpha=0.7)
plt.show()


#SELECT BACKGROUND GALAXIES


print '------------------------------------------------------------------------'
print '                     SELECTING BACKGROUND GALAXIES                      '
print '========================================================================'


#--------------------SELECT BACKGROUND GALAXIES------------------------------

#~ MAGMAX=moda+0.5
#~ MAGMIN=18.8
#~ print 'selecting galaxies between mag r ',MAGMIN,' and ',MAGMAX



sigmae=(erre1**2.0+erre2**2.0)**0.5


mask=np.full(ngxcat, True, dtype=bool)


#~ mask=(MAGBEST > MAGMIN)*(COLOR < 70.0)*(MAGBEST < MAGMAX)*(fwhm > 5.0)*(sigmae < 0.2)

faintgx=gxcat[mask,:]
faintim2=gxim2[mask,:]
faintmag=MAGBEST[mask]
faintcolor=COLOR[mask]
e1faint=e1[mask]
e2faint=e2[mask]
nf=len(faintcolor)
xfaint=posx[mask]				
yfaint=posy[mask]			
alfafaint=faintgx[:,3]
deltafaint=faintgx[:,4]

ro_gx=nf/(x.max()*y.max()*(pixsize**2)*((1.0/60.0)**2))
print 'densidad de gx = ', ro_gx, ' gx/arcmin**2'


		

print 'galaxias background ',nf

plt.plot(MAGBEST, COLOR, 'k.', faintmag, faintcolor, 'b.')
#~ plt.axis([15,30,-5,10])
plt.xlabel('i')
plt.ylabel('v-i')
plt.show()








idfaintsex=faintgx[:,0]
alfafaint=faintgx[:,3]
deltafaint=faintgx[:,4]
faint=np.zeros((nf,14),float)

# Asign a weigh to each galaxy given the probability that it's a background galaxy

# Asign a weigh to each galaxy given the probability that it's a background galaxy

prob_cat = np.loadtxt('prob_back_z'+str(z_cluster)+'.cat', comments='#')
peso=np.zeros(nf,float)

nbin=len(prob_cat)/5

print 'nbin', nbin

for i in range(5):
	colormin=prob_cat[i+nbin*i,0]
	colormax=prob_cat[i+nbin*i,1]	
	for k in range(nbin):
		magmin=prob_cat[k+nbin*i,2]
		magmax=prob_cat[k+nbin*i,3]
		for j in range(nf):
			if faintmag[j] < magmax and faintmag[j] > magmin and faintcolor[j] < colormax and faintcolor[j] >colormin:
				peso[j]=prob_cat[k+nbin*i,4]

plt.hist(peso, 40, facecolor='b', normed=1)
plt.title('distribucion de pesos')
plt.xlabel('peso')
plt.show()
#~ 
#~ plt.plot(e1faint,e2faint,'k.')
#~ plt.xlabel('e1')
#~ plt.ylabel('e2')
#~ plt.show()
#~ 


print '... and finally...'
print '------------------------------------------------------------------------'
print '             MAKING THE CATALOGUE FOR BACKGROUND GALAXIES               '
print '========================================================================'

faint[:,0]=idim2[mask]
faint[:,1]=faintgx[:,0]
faint[:,2]=posx[mask]
faint[:,3]=posy[mask]
faint[:,4]=faintgx[:,3]
faint[:,5]=faintgx[:,4]
faint[:,6]=faintmag
faint[:,7]=faintcolor
faint[:,8]=e1[mask]
faint[:,9]=e2[mask]
faint[:,10]=peso
faint[:,11]=sigmae[mask]
faint[:,12]=fwhm[mask]
faint[:,13]=ab[mask]
os.system('rm gx_back.cat')
f1=open('gx_back.cat','w')

f1.write('# Numero de galaxias background ')
f1.write(str(nf))
f1.write('\n')
f1.write('#1 SExtractor id \n')
f1.write('#2 Im2shape pos x \n')
f1.write('#3 Im2shape pos y \n')
f1.write('#5 Alfa J2000\n')
f1.write('#6 Delta J2000 \n')
f1.write('#7 MAGBEST \n')
f1.write('#8 COLOR \n')
f1.write('#9 e1 \n')
f1.write('#10 e2 \n')
f1.write('#11 peso \n')
f1.write('#12 sigma_e \n')
f1.write('#13 fwhm \n')
f1.write('#14 ab \n')
np.savetxt(f1, faint, fmt=['%6i']*2+['%12.9f']*12)
f1.close()



#LENSent file
#~ 
#~ lensent=np.zeros(((len(e1)),6),float)
#~ 
#~ e=(e1**2+e2**2)**0.5
#~ theta=(np.arctan2(e2,e1)/2.0)+(np.pi/2.0)
#~ 
#~ 
#~ 
#~ lensent[:,0]=(7.6419046-ra)*3600.
#~ lensent[:,1]=(26.302698-dec)*3600.
#~ lensent[:,2]=e*cos(2.*theta)
#~ lensent[:,3]=0.3
#~ lensent[:,4]=e*sin(2.*theta)
#~ lensent[:,5]=0.3
#~ 
#~ 
#~ 
#~ np.savetxt('input_lensent.cat', lensent, fmt=['%10.2f']*2+['%10.3f']*4)


print 'END OF PROGRAM :)'

