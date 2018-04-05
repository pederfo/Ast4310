from numpy import *
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sherif'})


def planck(T, wav):
	h = 6.62607E-27 #Planck constant
	c = 2.99792E10
	k = 1.380658E-16
	b = ((2.*h*c**2)/wav**5)*1/(exp( h*c/(wav*k*T) )-1)
	return b

wav = arange(1000,20801,200)
b = zeros(wav.shape)

plt.xlabel(r'wavelength $\lambda / \AA$', size=14)
plt.ylabel(r'Planck function', size=14)
plt.xlim(0,20800)
#plt.yscale('log')
plt.xscale('log')

for T in range(8000,5000-1,-200):
	b[:] = planck(T, wav[:]*1E-8)
	plt.plot(wav,b,'-')

plt.show()
