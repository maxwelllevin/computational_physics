# -*- coding: utf-8 -*-
"""
Exercise 6.9
With math to improve efficiency

Created on Thu Oct 19 17:50:16 2017

@author: Maxwell
"""

import numpy as np

pi = np.pi              # Pi
m_e = 9.1094e-31        # Mass of an electron
h = 6.6e-34             # Constant h
h_bar = h / (2 * pi)    # Constant 'h bar'
e = 1.6e-19             # Charge of an electron (1eV)
a = 1.6e-18             # Constant in our potential function (10eV)
L = 5e-10               # Width of our quantum well 
size = 100              # Size of H_mn array - marginal impact on accuracy


def H_mn(m,n):
    """ We do some math off-camera to simplify our large integral
        into a small piecewise function"""
    if m == n: 
        return (h_bar*h_bar*pi*pi*n*n) / (2*m_e*L*L) + a/2
        #(0.5/m_e)*((h_bar*pi*n/L)**2) + L*L/4
    elif (m%2 == n%2):
        return 0
    numer = -8*m*n*a
    denom = pi*pi*((m*m-n*n)**2)
    return numer/denom


def construct_H(size):
    """ Builds a square matrix H of a given size. """
    H = np.zeros([size,size], float)
    for m in range(size):
        for n in range(size):
            H[m,n] = H_mn(m+1,n+1)/e
    return H


def psi_prob_dist(x):
    total=0
    for i in range(size):
        total+=EIG[i]*np.sin(pi*(i+1)*x/L)
    return abs(total**2)



# Function that computes the integral of a given function from set bounds with specific depth (~accuracy)
def simpson_integral(lower_bound=0, upper_bound=L, function=lambda x: 1, depth=1e-14):
    N = int( (upper_bound - lower_bound) / depth )
    ans = ( function(lower_bound) + function(upper_bound))
    for k in range(1, N, 2):
        ans += 4 * function(lower_bound + k * depth)
    for k in range(2, N-1, 2):
        ans += 2 * function(lower_bound + k * depth)
    ans *= depth / 3
    return ans

# ================================== #

H = construct_H(size)
EIG, VEC = np.linalg.eigh(H)
distrib = simpson_integral(0, L, psi_prob_dist)

print("Our matrix H is:\n", H)
print("\nOur corresponding eigenvalues are:\n", EIG[:10])
print("Psi =", distrib)



