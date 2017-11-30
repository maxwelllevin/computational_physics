# -*- coding: utf-8 -*-
"""
HW #7 Problem 2

Created on Wed Nov 29 16:39:51 2017

@author: Maxwell Levin, with help from Iddo Fuhrmann
"""

import numpy as np
import pylab as pl

# Constants
M = 100         # Grid squares on a side
V = 1.0         # Voltage at top wall
delta = 1e-1    # Target accuracy
e0 = 1          # Permitivity of free space
r1 = 75
r2 = 25


# Create arrays to hold potential values
phi = np.zeros([M+1,M+1],float)
phiprime = np.zeros([M+1,M+1], float)


def eta(i,j):
    if i == 50: return 0
    return np.arctan((j-r2)/(i-50)) - np.arctan((j-r1)/(i-50))


# Main loop
error = 1.0
while error > delta:

    # Calculate new potential values
    for i in range(M+1):
        for j in range(M+1):
            if i == 0 or i == 50 or i == M or j == 0 or j == M:
                phiprime[i,j] = phi[i,j]
            else:
                phiprime[i,j] = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1] + np.sin(2*phi[i,j] + 2*eta(i,j)) ) / 4

    # Calculate Error
    error = np.max(np.abs(phi-phiprime))

    phi,phiprime = phiprime,phi

sigma = phi

pl.imshow(sigma)
#pl.show()



def dx(i, j):
    if (i != M-1):
        return sigma[i+1, j] - sigma[i,j]
    return sigma[i, j] - sigma[i-1, j]

def dy(i,j):
    if (j != M-1):
        return sigma[i,j+1] - sigma[i, j]
    return sigma[i, j] - sigma[i, j-1]


L = np.zeros([M,M], float)
for i in range(M):
    for j in range(M):
        L[i,j] = 0.5*( (dx(i,j)**2) + (dy(i,j)**2) + np.cos(2*sigma[i,j] + 2*eta(i,j)))

#pl.imshow(L)


