import sys
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from cosmolopy import *
from profiles_fit import *
cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.7}


gxcat1 = np.loadtxt('gx_back.cat', comments='#') #sextractor catalogue for galaxies

# EDIT THIS PARAMETERS
name_cluster='[VMF98]001'
filtro='r'
pixsize=0.1454
zcluster=0.501# pix size in arcsec
RIN=91.#bcg 90. e 100.    #where the fit starts, in kpc
ROUT=800.0 #fit limit in kpc
#~ BIN=10.0
beta=0.41
#~ err_beta=0.007
#~ ALFA0=168.3356    #BCG
#~ DELTA0=17.594405  #BCG 
#~ ALFA0=168.35984    #Estructura norte
#~ DELTA0=17.688599  #Estructura norte
ALFA0=168.537    #Estructura este
DELTA0=17.528045  #Estructura este

#~ x0=934  #BCG
#~ y0=908  #BCG
x0=935.#Estructura norte
y0=900.  #Estructura norte
#~ x0=1965  #Estructura este
#~ y0=2174  #Estructura este



MAGMIN=23.
MAGMAX=26.6
FWHM=5.0
SIGMAE=0.2

xy0=[x0,y0]

print 'Remember, if you did not edit the parameters in the program, this will not work!'

faintmag=gxcat1[:,6]-0.113
color=gxcat1[:,7]
fwhm=gxcat1[:,12]
sigmae=gxcat1[:,11]

ab=gxcat1[:,13]
e2=gxcat1[:,9]


mask=(faintmag > MAGMIN)*(color < 70.0)*(faintmag < MAGMAX)*(fwhm > 6.)*(sigmae < 0.2)*(ab<200.)

gxcat=gxcat1[mask,:]

#parameters
D_ang=cd.angular_diameter_distance(zcluster, z0=0, **cosmo)

kpcscale=D_ang*(((1.0/3600.0)*np.pi)/180.0)*1000.0

print 'kpcscale',kpcscale

cvel=299792458;   # Speed of light (m.s-1)
G= 6.67384e-11;   # Gravitational constant (m3.kg-1.s-2)
pc= 3.085678e16; # 1 pc (m)
Msun=1.989e30 # Solar mass (kg)


H=cd.hubble_z(zcluster,**cosmo) #H at z_cluster s-1
roc=(3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_cluster (kg.m-3)
Dl=D_ang*1.0e6*pc
#sigmac=(cvel**2.0)/(4.0*np.pi*G*Dl*beta)
sigmac=((cvel**2.0)/(4.0*np.pi*G*Dl))*(1./beta) #sup critical density (kg.m-2)
#~ err_sigmac=((cvel**2.0)/(4.0*np.pi*G*Dl))*(-1./(beta**2))*err_beta 
roc_mpc=roc*((pc*1.0e6)**3.0)
print 'densidad critica',(roc/((100.)**3))*1000.


Dcluster_m=D_ang*1.0e6*pc
sigmac=(cvel**2.0)/(4.0*np.pi*G*Dcluster_m*beta) #sup critical density (kg.m-2)
print 'sigmac',sigmac*(pc**2/Msun)

print 'n de gx', len(color)
posx=gxcat[:,2]
posy=gxcat[:,3]
e1=gxcat[:,8]
e2=gxcat[:,9]
sigmae=gxcat[:,11]
peso2=gxcat[:,10]
ra=gxcat[:,4]
dec=gxcat[:,5]
ab=gxcat[:,13]
fwhm=gxcat[:,12]
e=(e1**2+e2**2)**0.5
a=(ab*((1.0+e)/(1.0-e)))**0.5
theta=(np.arctan2(e2,e1)/2.0)+(np.pi/2.0)

#~ peso2=1./((fwhm.max()/fwhm)**2)
#~ peso2=peso2+peso
ngxcat=len(e)
area=posx.max()*posy.max()*((pixsize/60.)**2)

print 'densidad de galaxias background', len(e1)/(area)

#~ plt.plot(e1,e2,'k.')
#~ plt.show()

maska=a<15.
xbin,ybin,ex,ey,ngx=shear_map(posx[maska],posy[maska],a[maska],theta[maska],40)
#~ quiver(xbin,ybin,ex,ey,headwidth=1,headlength=0)
#~ plt.plot(x0,y0,'ro',x1,y1,'ro')

#~ xbin,ybin,ex,ey,ngx=shear_map2(posx,posy,e1,e2,7)
#~ ebin=np.sqrt(ex**2+ey**2)
#~ quiver(xbin,ybin,ex,ey,mode,cmap=cm.winter,headwidth=1,headlength=0,scale_units='xy',scale=0.0006)
#~ plt.colorbar()
#~ plt.plot(x0,y0,'bo',xe,ye,'bo')
#~ plt.show()

#~ thetabin=(np.arctan2(ex,ey)/2.0)+(np.pi/2.0)
#~ quiver(xbin,ybin,np.cos(2.*thetabin),np.sin(2.*thetabin),ebin,cmap=cm.winter,headwidth=1,headlength=0,scale_units='xy',scale=0.0005)
#~ plt.colorbar()
#~ plt.plot(x0,y0,'bo',xe,ye,'bo')
#~ plt.show()
#~ quiver(xbin,ybin,ebin*np.cos(2.*thetabin),ebin*np.sin(2.*thetabin),ebin,cmap=cm.winter,headwidth=1,headlength=0,scale_units='xy',scale=0.01)
#~ plt.colorbar()
#~ plt.plot(x0,y0,'ro')
#~ plt.show()


def shear_profile_log(RIN,ROUT,r,et,ex,peso,STEP):
	
	nbin=15
	
	SHEAR=np.zeros(nbin,float)
	CERO=np.zeros(nbin,float)
	R=np.zeros(nbin,float)
	err=np.zeros(nbin,float)
	
	stepbin=STEP
	rin=RIN/kpcscale#120
	rout=10**(np.log10(rin)+stepbin)
	BIN=0
	
	for j in range(nbin):
		if (rin+rout)/2.0 < ROUT/kpcscale:
			BIN+=1
		maskr=(r>rin)*(r<rout)	
		w2=peso[maskr]
		pes2=w2.sum()			
		shear=et[maskr]*w2
		cero=ex[maskr]*w2
		ERR=w2**2
		n=len(shear)
		R[j]=rin+(rout-rin)/2.0	
		if n == 0:
			SHEAR[j]=0.0
			CERO[j]=0.0
			err[j]=0.0
		else:	
			SHEAR[j]=(shear.sum()/pes2)
			CERO[j]=cero.sum()/pes2
			#~ err[j]=(0.3/(n**0.5))
			sigma_e=(0.3**2.)#+dife[maskr]**2.)
			ERR2=ERR.sum()
			err[j]=((ERR2*sigma_e)**0.5)/pes2
		#print 'shear',n, rin, rout, SHEAR[j], err[j]
		rin=rout
		rout=10**(np.log10(rin)+stepbin)

	#~ plt.hist(dife,30)
	#~ plt.xlabel('DIFE')
	#~ plt.show()
	#print 'galaxias usadas', contar, BIN
	
	
	return [R,SHEAR,CERO,err,BIN]

def radianes(x):
	y=(x*np.pi)/180.
	return y
	
def grados(x):
	y=(x*180.)/np.pi
	return y

def rot(posx,posy,e1,e2,xy0):

		#print '------------------------------------------------------------------------'
		#print '                     CHANGING THE CORDINATE AXIS                        '
		#print '========================================================================'
		ngxcat=len(posx)
		et=np.zeros(ngxcat,float)
		ex=np.zeros(ngxcat,float)	
		e=(e1**2+e2**2)**0.5
		theta=(np.arctan2(e2,e1)/2.0)+(np.pi/2.0)
		print 'ngxcat',ngxcat, len(posx)
		for j in range(ngxcat):
			u1=posx[j]-xy0[0]
			u2=posy[j]-xy0[1]
			if u1 != 0.0 and u2 != 0.0:
				fi=np.arctan2(u2,u1)
				if u2 < 0.0:    #check if the angle is in III or IV quadrant
					fi=2.0*np.pi+fi
				beta2=theta[j]-fi#+(np.pi/4.0)
				et[j]=(-1.0*e[j]*cos(2.0*beta2))
				ex[j]=(e[j]*sin(2.0*beta2))
				
		return et,ex


#~ r=(grados(np.arccos(np.sin(radianes(dec))*np.sin(radianes(DELTA0))+np.cos(radianes(dec))*np.cos(radianes(DELTA0))*np.cos(radianes(ra-ALFA0))))*3600.)

r=((((posx-xy0[0])**2+(posy-xy0[1])**2)**0.5)*pixsize)


et,ex=rot(posx,posy,e1,e2,xy0)
#~ peso2.fill(1.)

print '---------------------------------------------------------'
print '      COMPUTING THE STEP ACORDING THE ERROR IN DISP      '
print '========================================================='

param=1000000000.0
for paso in range(30):

	stepbin=0.14+0.001*paso #bcg 0.14
	print 'paso',stepbin
	
	R,shear,CERO,err,BIN=shear_profile_log(RIN,ROUT,r,et,ex,peso2,stepbin)
	#compute the Einstein ring for a SIS model
	R=(R*kpcscale)/1.e3 #R en Mpc
	sis=SIS_fit(R[:BIN],shear[:BIN],err[:BIN],beta,zcluster)
	#~ sis=NFW_fit(R[:BIN],shear[:BIN],err[:BIN],4.,zcluster,sigmac)
	disp=sis[0]
	errordisp=sis[1]
	chired=sis[2]
	sigma=disp/errordisp
	print 'error en sigma',sigma,' disp',disp,'+/-',errordisp, BIN
	
	if chired < param:
		param=chired
		STEP=stepbin

try:
	print 'STEP',STEP
except:
	STEP=paso
# EDIT THIS PART TOO







print '---------------------------------------------------------'
print '             COMPUTING THE SHEAR PROFILES                '
print '========================================================='


R,SHEAR,CERO,err,BIN=shear_profile_log(RIN,ROUT,r,et,ex,peso2,STEP)
R=(R*kpcscale)/1.e3 #R en Mpc
print '---------------------------------------------------------'
print '                   FITTING PROFILES                      '
print '========================================================='


print 'First a SIS profile'



shear=SHEAR
cero=CERO

sis=SIS_fit(R[:BIN],shear[:BIN],err[:BIN],beta,zcluster)

rin=0.0
rout=ROUT/(kpcscale)
lim=rout
print 'limite ajuste',rout


disp=sis[0]
errordisp=sis[1]



CHI_sis=sis[2]





x1=sis[3]
y1=sis[4]



X=x1
Y=y1
DISP=disp
ERRORDISP=errordisp

x=np.arange(0,10000)
y=np.zeros(len(x))



rin=(rin*kpcscale)/1000.0 
rout=(rout*kpcscale)/1000.0 




# M200 SIS

M200_SIS=((2.*(disp*1.e3)**3)/((50**0.5)*G*H))/(Msun)
e_m200_SIS=(((6.*(disp*1.e3)**2)/((50**0.5)*G*H))*(errordisp*1.e3))/(Msun)

print 'M200(SIS) =', '%.1e' % M200_SIS, '+/-','%.1e' % e_m200_SIS



print 'Now is trying to fit a NFW profile...'
print 'fitting c...'

        
c=4.
    
    
#try to fit NFW profile (ecuaciones del paper de King and Schneider 2001)

print 'c =',int(c)

nfw=NFW_fit(R[:BIN],shear[:BIN],err[:BIN],c,zcluster,sigmac)

CHI_nfw=nfw[2]
print 'chi cuadrado reducida NFW = ', CHI_nfw


R200=nfw[0]

error_R200=nfw[1]

print 'R200 =',R200,'+/-', '%.2f' % error_R200,'Mpc'

#make the plot
x2=nfw[3]
y2=nfw[4]
#Compute the mass
# M200 mass

M200_NFW=(800.0*np.pi*roc_mpc*(R200**3))/(3.0*Msun)

e_M200_NFW=((800.0*np.pi*roc_mpc*(R200**2))/(Msun))*error_R200
print 'M200 nfw =',M200_NFW,'+/-',e_M200_NFW
MR200_SIS=(2.*((DISP*1.e3)**2)*(R200*1.0e6*pc))/(G*Msun)
e_MR200_SIS=((4.*(DISP*1.e3)*(R200*1.0e6*pc))/(G*Msun))*(ERRORDISP*1e3)

print 'M200 sis =',MR200_SIS,'+/-',e_MR200_SIS


if ERRORDISP==inf  or error_R200==inf:
	ERRORDISP=999.
	error_R200=999.

f2=open('../../result_ajuste','a')
tabla=np.zeros((1,11),float)
tabla[0,0]=001
tabla[0,1]=CHI_sis
tabla[0,2]=round(DISP,-1)
tabla[0,3]=round(ERRORDISP,-1)
tabla[0,4]=(M200_SIS)/1.e14
tabla[0,5]=(e_m200_SIS)/1.e14
tabla[0,6]=CHI_nfw
tabla[0,7]=R200
tabla[0,8]=error_R200
tabla[0,9]=(M200_NFW)/1.e14
tabla[0,10]=(e_M200_NFW)/1.e14

np.savetxt(f2, tabla, fmt=['%4i']*1+['%6.1f']*1+['%4i']*2+['%6.1f']*7)
f2.close()


# ------- FOR THE PLOT ------

x=np.zeros(2,float)
y=np.zeros(2,float)
X2=X
x[1]=6

blancox=1000.
blancoy=1000.0

rcParams['font.family'] = 'serif'
rcParams['figure.figsize'] = 9.5, 3.0
fig, ax = plt.subplots()
majorFormatter = FormatStrFormatter('%.1f')

plt.plot(X,Y,'k',label='hola')
plt.axis([RIN/kpcscale,lim,-0.4,1.0])
plt.legend()
matplotlib.rcParams['legend.fontsize'] = 12.0
#~ plt.xlabel(u'r [arcsec]',fontsize=13)
plt.ylabel(r'$\langle \gamma \rangle$')
#~ plt.suptitle(name_cluster,fontsize=13.5)
plt.twiny()
ax.yaxis.set_major_formatter(majorFormatter)
plt.plot(blancox,blancoy,'w.',label=name_cluster)
plt.plot(R,shear,'ko')
plt.plot(R,cero,'kx')
plt.plot(x2,y2,'k--',label='NFW: $R_{200}$ = '+str('%.1f' % R200)+' $\pm$ '+str('%.1f' % error_R200)+' Mpc$\,h^{-1}_{70}$, $\chi_{red}^{2} =$'+str('%.1f' % CHI_nfw)) #, c='+str(round(c)))
plt.plot(x,y,'k')
plt.plot(X2,Y,'k',label='SIS: $\sigma$ = '+str(int(round(DISP,-1)))+' $\pm$ '+str(int(round(ERRORDISP,-1)))+' km/s, $\chi_{red}^{2} =$'+str('%.1f' % CHI_sis))
matplotlib.rcParams['legend.fontsize'] = 10.0
plt.legend()
plt.errorbar(R, shear, yerr=err, fmt=None, ecolor='k')
plt.errorbar(R,cero, yerr=err, fmt=None, ecolor='k')
plt.xlabel('r [Mpc$\,h^{-1}_{70}$]')
plt.ylabel('Shear')
plt.axis([RIN*1.0e-3,(lim*kpcscale*1.0e-3),-0.1,0.4])
ax.yaxis.set_ticks(np.arange(-0.1, 0.4, 0.1))
plt.subplots_adjust(left=0.1, right=0.85, top=0.85, bottom=0.1)
plt.savefig('shear_profile_'+name_cluster+'.eps', format='eps',bbox_inches='tight')
plt.show()


print 'END OF PROGRAM :)'
