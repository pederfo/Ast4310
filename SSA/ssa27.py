from numpy import *
import matplotlib.pyplot as plt


#Populations of hydrogen

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
	nstage[1] = nstage[0] * sahaconst * u[1]/u[0] * exp( -13.598/KeV )
	ntotal = sum(nstage)



	#Boltzmann

	nlevel = nstage[0]*g[0,level-1]/u[0]*exp( -chiexc[0,level-1]/kevT )
	nlevelrel = nlevel/ntotal

	for s in range(6):
		print s+1, g[0,s], chiexc[0,s], g[0,s]*exp( -chiexc[0,s]/kevT )

	for s in range(0, nrlevels, 10):
		print s+1, g[0,s], chiexc[0,s], g[0,s]*exp( -chiexc[0,s]/kevT )

	return nlevelrel

print sahabolt_H(6000,1.E2,1.)




	
	
