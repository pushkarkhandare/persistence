import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

delta = 0.2
eps = 1.0
sigma = 1.0
a = 1.0
width = 20 
height = 19
N = (width*(height + 1)/2) + ((width-1)*(height - 1)/2)
lat = np.zeros((N,2))

for k in range(0,height):
    if k%2==0:
        for i in range((k*width) - (k/2) ,((k+1)*width) - (k/2)):
            lat[i] = [a*(i-((k*width)-(k/2))),np.sqrt(3)*k/2]
    if k%2==1:
        for i in range((k*width) - ((k-1)/2),(k+1)*(width - (1/2))):
            lat[i] = [(a/2) + a*(i-((k*width)-((k-1)/2))),k*np.sqrt(0.75)*a]

def distfunc(x0,x1):
    dimensions = np.array([width,height])
    delta = np.abs(x0-x1)
    delta = np.where(delta>0.5*dimensions,delta-dimensions,delta)
    return np.sqrt((delta**2).sum(axis=-1))

print distfunc([[1,2],[2,3,[4,6]],[[0,3],[2,3],[4,4]])

def ljp(r):
        return 4*(np.exp(-12*np.log(r)) - np.exp(-6*np.log(r)))

#def mcmove(delta):
#    rand = delta*np.random.random(2)
#    rand = np.where(rand>delta*0.5,rand-delta, rand)
#    return rand
 

#def neigh(k):
#    Neigh = []
#    for i in range(0,N):
#        dist = distfunc(lat[k,0],lat[k,1],lat[i,0],lat[i,1])
#        if dist<3.0 and dist!=0:
#            Neigh.append([lat[i,0],lat[i,1]])
#    return np.array(Neigh)
#
#def neighlist():
#    neighlist = []
#    for i in range(0,N):
#        neighlist.append(neigh(i))
#    return np.array(neighlist)

#def energy_calc(k):
#    A = neighlist()[k]
#    n = len(A)
#    energy = 0
#    for i in range(0,n):
#        dist = distfunc(lat[k,0],lat[k,1],A[i][0],A[i][1])
#        energy = energy + ljp(dist)
#    return energy

def direct_energy(k):
    dists = distfunc(lat[k],lat)
    relevant = ma.masked_where(dists>3,dists)
    relevant = ma.masked_where(relevant==0, relevant)
    energies = ljp(relevant)
    return energies.sum()

accept=0
#reject=0
k=0
mcenergy = np.zeros(5000)
#while k<5000:
#    print k
#    for i in range(0,N):
#        rand_index = np.random.randint(0,N)
#        Eold = direct_energy(rand_index)
#        rand = delta*np.random.random(2)
#        rand = np.where(rand>delta*0.5,rand-delta, rand)
#        lat[rand_index,0] = lat[rand_index,0]+rand[0]
#        lat[rand_index,1] = lat[rand_index,1]+rand[1]
#        Enew = direct_energy(rand_index)
#        Ediff = Enew - Eold
#        if Ediff<0: 
#            accept = accept + 1
#            mcenergy[k] = mcenergy[k] + Enew
#        if Ediff>0:
#            p = np.random.random()
#            if p < np.exp(-Ediff):
#                accept = accept + 1
#                mcenergy[k] = mcenergy[k] + Enew
#            else:
#               # reject = reject + 1
#                lat[rand_index,0] = lat[rand_index,0]-rand[0]
#                lat[rand_index,1] = lat[rand_index,1]-rand[1]
#                mcenergy[k] = mcenergy[k] + Eold
#    k = k+1
#
#print accept/(5000.0*N)
#plt.plot(mcenergy)
#plt.show()
