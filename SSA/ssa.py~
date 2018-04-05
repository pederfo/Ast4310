from numpy import *
import matplotlib.pyplot as plt

#constants
#temp = 20000.

def u_func(temp):
	chiion = array([7, 16, 31, 51]) #Schadee ionization energies
	k = 8.61734E-5 #Boltzmann constant
	u = zeros(4) #partition function 
	for r in range(len(u)):  #Fill u vector; U1, U2 = U3 = U4
		for s in range(chiion[r]):
			u[r] = u[r] + exp( -s/k/temp )
	return u

def boltz(temp, r, s):
	u = u_func(temp)
	KeV = 8.61734E-5
	relnrs = 1./u[r-1]*exp(-(s-1)/(KeV*temp))
	return relnrs

def saha(temp, elpress, ionstage):
	kerg = 1.380658E-16
	KeV = 8.61734E-5 # Boltzmann constant
	h = 6.62607E-27 #Planck constant
	elmass = 9.109390E-28 # electron mass
	kevT = KeV * temp
	kergT = kerg * temp
	eldens = elpress / kergT
	chiion = array([7, 16, 31, 51])
	u = u_func(temp)
	u = append(u, 2)
	sahaconst = (2.*pi*elmass*kergT/(h**2))**1.5*2./eldens
	nstage = zeros(5)
	nstage[0] = 1.
	
	for r in range(4):
		nstage[r + 1] = nstage[r]*sahaconst*u[r +1]/u[r]*exp(-chiion[r]/kevT)
	ntotal = sum(nstage)
	nstagerel = nstage/ntotal

	return nstagerel[ionstage -1]


def sahaboltz(temp, elpress, ion, level):
	
	return saha(temp, elpress, ion) * boltz(temp, ion, level)




temp = arange(0, 30001, 1000)
pop = zeros((5, 31))
for T in arange(1,31):
	for r in arange(1,5):
		pop[r,T] = sahaboltz(temp[T], 131., r, 1.)

plt.figure(0)
for i in range(1,5):
	plt.plot(temp,pop[i,:], label=i)
plt.title('r = 1, s = 1')
plt.xlabel('temperature', size=14)
plt.ylabel('population', size=14)
plt.yscale('log')
plt.ylim([1E-3,1.1])
plt.legend(loc=4, fontsize=12)
plt.show()

temp = arange(0, 30001, 1000)
pop = zeros((5, 31))
for T in arange(1,31):
	for r in arange(1,5):
		pop[r,T] = sahaboltz(temp[T], 131., r, 2.)

plt.figure(1)
for i in range(1,5):
	plt.plot(temp,pop[i,:], label=i)
plt.title('r = 1, s = 2')
plt.xlabel('temperature', size=14)
plt.ylabel('population', size=14)
plt.yscale('log')
plt.ylim([1E-3,1.1])
plt.legend(loc=4, fontsize=12)
plt.show()

temp = arange(0, 30001, 1000)
pop = zeros((5, 31))
for T in arange(1,31):
	for r in arange(1,5):
		pop[r,T] = sahaboltz(temp[T], 131., r, 4.)

plt.figure(2)
for i in range(1,5):
	plt.plot(temp,pop[i,:], label=i)
plt.title('r = 1, s = 4')
plt.xlabel('temperature', size=14)
plt.ylabel('population', size=14)
plt.yscale('log')
plt.ylim([1E-3,1.1])
plt.legend(loc=4, fontsize=12)
plt.show()


