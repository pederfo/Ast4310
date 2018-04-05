from numpy import *
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sherif'})
#from scipy import special
import math


def voigt(a,u):
	lamd = 1.
	return 1./(lamd*sqrt(pi))*(exp(-u**2)+ a/(sqrt(pi)*u**2))

def planck(T, wav):
	h = 6.62607E-27 #Planck constant
	c = 2.99792E10
	k = 1.380658E-16
	b = ((2.*h*c**2)/wav**5)*1/(exp( h*c/(wav*k*T) )-1)
	return b



u = arange(-10.,10.1,0.1)
a = array([0.001, 0.01, 0.01, 1.])
vau = zeros((a.shape[0], u.shape[0]))

for i in range(4):
	vau[i,:] = voigt(a[i],u[:])
	plt.plot(u[:], vau[i,:], label = 'a = ' +str(a[i]))

#plt.ylim(0,1)
plt.xlim(-10,10)
plt.legend(fontsize=12)
plt.ylabel('voigt profile', size=12)
plt.yscale('log')
#plt.xscale('log')


## Schuster-Schwarszchild line profile

Ts = 5700.
#Ts = 10200.
T1 = 4200.
#T1 = 3200.
a = 0.1
wav = 5000.0e-8
tau0 = 1.
u = arange(-10.,10.1,0.1)
intensity = zeros(u.shape)

for i in range(201):
	tau = tau0*voigt(a,u[i])
	intensity[i] = planck(Ts,wav)*exp(-tau) + planck(T1,wav)*(1.-exp(-tau))

plt.plot(u,intensity)
plt.show()


logtau0 = arange(-20.,20.1,0.5)

for itau in range(9):
	for i in range(201):
		tau = 10.**(logtau0[itau] * voigt(a,u[i]))
		intensity[i] = planck(Ts,wav)*exp(-tau) + planck(T1,wav)*(1.-exp(-tau))
	plt.plot(u,intensity, label = r'$\log{(\tau_0)} = $' + str(logtau0[itau]))

#plt.legend(loc=3, fontsize=12)
plt.ylabel('Intensity')
plt.xlabel('Wavelength (u)')
plt.show()



for iwav in range(1,4):
	wav = (iwav*2+1.)*1.0e-5
	for itau in range(8):
		for i in range(201):
			tau = 10.**(logtau0[itau]) * voigt(a,u[i])
			intensity[i] = planck(Ts,wav)*exp(-tau) + planck(T1,wav)*(1.-exp(-tau))
		intensity = intensity / intensity[0]
		plt.plot(u,intensity[:], linewidth=1.)


plt.ylabel('Intensity')
plt.xlabel('Wavelength (u)')
plt.show()
