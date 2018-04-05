from numpy import *
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sherif'})

falc = loadtxt("falc.dat")
#print lines
#print len(lines)

#Various arrays from falc.dat file

h = falc[0,:]
tau5 = falc[1,:]
colm = falc[2,:]
temp = falc[3,:]
vturb = falc[4,:]
nhyd = falc[5,:]
nprot = falc[6,:]
nel = falc[7,:]
ptot = falc[8,:]
pgasptot = falc[9,:]
dens = falc[10,:]

#plot falc model
"""
plt.plot(h, temp)
plt.xlabel('height [km]')
plt.ylabel('temperature [k]')
"""
plt.plot(ptot, colm)
plt.yscale('log')
plt.xscale('log')



plt.show()




