from numpy import *
import matplotlib.pyplot as plt

#constants
#temp = 20000.

def u_func(temp):
	chiion = array([6, 12, 51, 67]) #Schadee ionization energies
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
	chiion = array([6, 12, 51, 67])
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

def sahabolt_H(temp, elpress, level):
	kerg = 1.380658E-16
	KeV = 8.61734E-5 # Boltzmann constant
	h = 6.62607E-27 #Planck constant
	elmass = 9.109390E-28 # electron mass
	kevT = KeV * temp
	kergT = kerg * temp
	eldens = elpress / kergT
	
	#energylevels and weights for hydrogen
	nrlevels = 100 			#partition function cutoff
	g = zeros((2,nrlevels))		#decleration weights
	chiexc = zeros((2,nrlevels))	#decleration excitation energies

	for s in range(nrlevels):
		g[0,s] = 2.*(s+1.)**2.			#statistical weights
		chiexc[0,s] = 13.598*(1.-1./(s+1.)**2.)	#excitation weights

	g[1,0] = 1.			#statistical weights free proton
	chiexc[1,0] = 0

	#partition functions

	u = zeros([2])
	for s in range(nrlevels):
		u[0] = u[0] + g[0,s]*exp( -chiexc[0,s]/KeV )
	u[1] = g[1,0]

	#saha

	sahaconst = (2.*pi*elmass*kergT/(h**2))**1.5*2./eldens
	nstage = zeros(2)
	nstage[0] = 1
	nstage[1] = nstage[0] * sahaconst * u[1]/u[0] * exp( -13.598/kevT )
	ntotal = sum(nstage)



	#Boltzmann
	
	nlevel = nstage[0]*g[0,level-1]/u[0]*exp( -chiexc[0,level-1]/kevT )
	nlevelrel = nlevel/ntotal
	"""
	for s in range(6):
		print s+1, g[0,s], chiexc[0,s], g[0,s]*exp( -chiexc[0,s]/kevT )

	for s in range(0, nrlevels, 10):
		print s+1, g[0,s], chiexc[0,s], g[0,s]*exp( -chiexc[0,s]/kevT )
	"""
	return nlevelrel



temp = arange(1000, 20001, 100)
CaH = zeros(temp.shape)
Caabund = 2.0E-6
for i in range(0,191):
	NCa = sahaboltz(temp[i],1E2, 2, 1)
	NH = sahabolt_H(temp[i],1E2, 2)
	CaH[i] = NCa*Caabund/NH

plt.plot(temp,CaH)
plt.title('strength ratio Ca II K / H alpha')
plt.yscale('log')
plt.xlabel('temperature T / K')
plt.ylabel('Ca II K / H alpha')
plt.show()


print 'CaH/H ratio at 5000 K = ', CaH[argwhere(temp==5000)][0][0]

## 2.9: Temperature sensitivity

temp = arange(2000,12001,100)
dNCadT = zeros(temp.shape)
dNHdt = zeros(temp.shape)
dT = 1.

for i in range(101):
	NCa = sahaboltz(temp[i],1E2,2,1)
	NCa2 = sahaboltz(temp[i]-dT,1E2,2,1)
	dNCadT[i] = (NCa - NCa2)/(dT*NCa)
	NH = sahabolt_H(temp[i],1E2,2)
	NH2 = sahabolt_H(temp[i]-dT,1E2,2)
	dNHdt[i] = (NH-NH2)/(dT*NH)

plt.figure()
plt.plot(temp,absolute(dNHdt), label='H')
plt.plot(temp,absolute(dNCadT), label='Ca II K')

plt.yscale('log')
plt.xlabel('temperature T / K')
plt.ylabel('(dn/dT)/n')

NCa = zeros(temp.shape)
NH = zeros(temp.shape)

for i in range(101):
	NCa[i] = sahaboltz(temp[i],1E2,2,1)
	NH[i] = sahabolt_H(temp[i],1E2,2)

plt.plot(temp,NH/amax(NH), ls='--', label = 'rel. pop. H')
plt.plot(temp,NCa/amax(NCa), ls='--', label = 'rel. pop. Ca II K')
plt.legend(loc=4, fontsize=12)
plt.show()


## 2.10 Hot stars versus cold stars

for T in arange(2E3,2E4+1.,2E3):
	print T, sahabolt_H(T,1E2,1.)

temp = arange(1E3, 2E4+1., 1E2)
nH = zeros(temp.shape)
for i in range(191):
	nH[i] = sahabolt_H(temp[i],1E2,1.)

plt.figure()
plt.plot(temp,nH)
plt.xlabel('temperature T / K', size=14)
plt.ylabel('neutral hydrogen fraction', size=14)
plt.legend()
plt.show()


