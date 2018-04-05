from numpy import *
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sherif'})


B = 2.
tau = arange(0.01,10.01,0.01)
#d = zeros(tau.shape)
intensity = zeros(tau.shape)
for IO in range(4,-1,-1):
	intensity[:] = IO * exp(-tau[:]) + B*(1-exp(tau[:]))
	#d[:] = 2
	#plt.plot(tau,  d)
	plt.plot(tau,intensity, label = 'intensity IO = ' +str(IO))


plt.xlabel(r'optical depth $ \tau$', size = 14)
plt.ylabel('intensity', size=14)
plt.legend(loc=4,fontsize=12)
plt.yscale('log')
plt.xscale('log')
plt.show()

