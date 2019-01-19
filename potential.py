import numpy as np
import matplotlib.pyplot as plt

r = np.arange(0.45,2,0.001)
pot = np.zeros(len(r))
sigma = 0.5
r0 = 1.0
d = 2.0
k = (2*d)/((sigma - r0)**2)

for i in range(0,len(r)):
    if r[i] < sigma:
        pot[i] = (sigma/r[i])**12 - (sigma/r[i])**6
    if r[i]>=sigma and r[i]<2*r0 - sigma:
        pot[i] = k*0.5*((r[i] - r0)**2) - d
    if r[i]>=2*r0 - sigma:
        pot[i] = 0
plt.xlabel('r',fontsize = 'xx-large')
plt.ylabel('Potential',fontsize='xx-large')
plt.plot(r,pot)
plt.axhline(y=0, color='r',linestyle = '--') 
plt.show()

