import numpy as np
import matplotlib.pyplot as plt

A=np.loadtxt('triang.dat')

a = A[0][0]
i=1
while A[i][0]!=a:
    i = i + 1
width = (i + 1)/2


