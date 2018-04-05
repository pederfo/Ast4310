from numpy import *
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sherif'})

def voigt(a,u):
	lamd = 1.
	return 1/(lamd*sqrt(pi))*(exp(-u**2) + a/(sqrt(pi)*u**2))

def planck(T, wav):
	h = 6.62607E-27 #Planck constant
	c = 2.99792E10
	k = 1.380658E-16
	b = ((2.*h*c**2)/wav**5)*1/(exp( h*c/(wav*k*T) )-1)
	return b




## 3.4 The equivalent with of spectral lines

def profile(a,tau0,u):
	Ts = 5700.
	T1 = 4200.
	wav = 5000.0e-8
	intensity = zeros(u.size)
	usize = u.size
	for i in range(usize):
		tau = tau0 * voigt(a,u[i])
		intensity[i] = planck(Ts,wav)*exp(-tau) + planck(T1,wav)*(1.-exp(-tau))
	return intensity

#Checking the profile
u = arange(-200,200.4,0.4)
a = 0.1
tau0 = 1.0e2
intensity = profile(a,tau0,u)

plt.plot(u,intensity)
plt.xlabel('u', size=14)
plt.ylabel('Intensity', size=12)

plt.show()

#relative
reldepth = (intensity[0]-intensity)/intensity[0]
plt.plot(u,reldepth)
plt.xlabel('u', size=14)
plt.ylabel(r'equivalent width $W_{\lambda}$', size=14)

eqw = sum(reldepth)*0.4
print eqw
plt.show()




## 3.5 The curve of growth

tau0 = logspace(-2,4,61)
eqw = zeros(tau0.size)

for i in range(61):
	intensity = profile(a,tau0[i],u)
	reldepth = (intensity[0] - intensity) / intensity[0]
	eqw[i] = sum(reldepth)*0.4

plt.plot(tau0,eqw)
plt.xlabel(r'$\tau_0$', size=18)
plt.ylabel(r'equivalent width $W_{\lambda}$', size=14)
plt.xscale('log')
plt.yscale('log')
plt.show()


