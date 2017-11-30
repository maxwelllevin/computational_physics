# -*- coding: utf-8 -*-
"""
Exercise 9.1: Poisson's Equation

Created on Thu Nov 23 16:39:40 2017

@author: Maxwell Levin
"""

import numpy as np
import pylab as pl

# Constants
M = 100         # Grid squares on a side
V = 1.0         # Voltage at top wall
delta = 1e-6    # Target accuracy
e0 = 1          # Permitivity of free space

# Create arrays to hold potential values
phi = np.zeros([M+1,M+1],float)
#phi[20:40,20:40] = -V
#phi[60:80,60:80] = V
phiprime = np.zeros([M+1,M+1], float)


def ro(i,j):
    if (i <= 80 and i >= 60 and j <= 40 and j >= 20):
        return 1
    elif (i <= 40 and i >= 20 and j <= 80 and j >= 60):
        return -1
    else: return 0

# Main loop
error = 1.0
while error > delta:

    # Calculate new potential values
    for i in range(M+1):
        for j in range(M+1):
            if i == 0 or i == M or j == 0 or j == M:
                phiprime[i,j] = phi[i,j]
            else:
                phiprime[i,j] = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1] + ro(i,j)) / 4

    # Calculate Error
    error = np.max(np.abs(phi-phiprime))

    phi,phiprime = phiprime,phi


pl.imshow(phi)
pl.show()