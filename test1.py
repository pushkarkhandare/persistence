import numpy as np
import matplotlib.pyplot as plt

A = np.loadtxt("energylow.dat")
A = A[500:]
tmax = len(A)
B = np.zeros(len(A))

A = A-np.mean(A)
Anorm = np.sum(A**2)
acor = np.correlate(A, A, "same")/Anorm
# use only second half
acor = acor[len(acor)/2:]
plt.xlabel("MCSteps",size='14')
plt.ylabel("Energy Autocorrelation",size='14')
plt.semilogx(acor)
plt.show()
#corr = np.correlate(A,A,mode='full')
#corr = corr[corr.size/2:]
#corr2 = np.mean(A)**2
#corr3 = A**2
#corr3 = np.mean(corr3)
#corr4 = corr3 - corr2
#print corr4
#corr4 = corr3 - corr2
#corr4i = 1.0/corr4
#corr = corr*corr4i - corr2*corr4i
#corr = [(i - corr2)/corr4 for i in corr]
#print corr[0]
#plt.plot(corr)
#plt.show()
