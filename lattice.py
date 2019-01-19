import numpy as np
import matplotlib.pyplot as plt

a = 0.5
L = 5.0
gamma = 0.0

rows = int(np.round(2*L*(3**-0.5)*(a**-1)))
H = rows*0.5*(3**0.5)*a
n_per_row = int(L/a)

N = rows*n_per_row
offset = np.array([(-L/2.0) + (a/2.0),(-H/2.0) + ((3**0.5)*0.25*a)])
print L,H,rows,n_per_row,N

lat = np.zeros((N,2))
for k in range(0,rows):
        for i in range(k*n_per_row,(k+1)*n_per_row):
            count = i - k*n_per_row
            if k%2==0:
                lat[i] = [a*count,k*0.5*(3**0.5)*a] + offset
            else:
                lat[i] = [a*(count - 0.5),k*0.5*(3**0.5)*a] + offset
                
if gamma!=0:
    print gamma
    for i in range(0,N):
        lat[i,0] = lat[i,0] + gamma*lat[i,1]


randlat = np.zeros((N,2))
#print H
for i in range(0,N):
    randlat[i,0] = (n_per_row-1.0)*(0.5 - np.random.rand())
    randlat[i,1] = H*(0.5 - np.random.rand())

#plt.scatter(randlat[:,0],randlat[:,1])
#plt.show()
def distfunc(x1,y1,x2,y2):
    return ((x2-x1)**2 + (y2 - y1)**2)**0.5

def neighlist(x,y):
    A=[]
    for i in range(0,N):
        if distfunc(x,y,lat[i,0],lat[i,1])<=1.2*a:
            A.append(lat[i])
    return A

outfile = open('latticedata.dat', 'w')
f = open('triang.dat','w')
f.write("%2.6f %2.6f %i %1.2f\n" % (L,H,N,a))
g = open('random.dat','w')

for i in range(0,N):
    f.write("%4.4f %4.4f\n" % (lat[i,0],lat[i,1]))
    g.write("%4.4f %4.4f\n" % (randlat[i,0],randlat[i,1]))
    B = neighlist(lat[i,0],lat[i,1])
    for j in range(0,len(B)):
        outfile.write("%4.4f %4.4f %4.4f %4.4f\n" % (lat[i,0],lat[i,1],B[j][0],B[j][1]))
outfile.close()
f.close()
g.close()
