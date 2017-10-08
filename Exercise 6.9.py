# -*- coding: utf-8 -*-
"""
Exercise 6.9

Created on Wed Oct  4 15:57:51 2017

@author: Maxwell Levin
"""

import numpy as np

# Declare our constants to use 
pi = np.pi
m_e = 9.1e-31 #kg
h = 6.6e-34 #J
h_bar = h / (2 * pi) #J
e = 1.6e-19 #J
a = 1.6e-18 #J
L = 5e-10 #m

# Returns the energy at a specific n-state
def find_E_n(n):
    return ( n*n * pi*pi * h_bar*h_bar ) / ( 2 * m_e * L*L )


# Function that computes the integral of a given function from set bounds with specific depth (~accuracy)
def simpson_integral(lower_bound=0, upper_bound=L, function=lambda x: 1, depth=0.1):
    N = int( (upper_bound - lower_bound) / depth )
    ans = ( function(lower_bound) + function(upper_bound))
    
    for k in range(1, N, 2):
        ans += 4 * function(lower_bound + k * depth)
    
    for k in range(2, N-1, 2):
        ans += 2 * function(lower_bound + k * depth)
    
    ans *= depth / 3
    return ans


# function inside the integral
def H_mn_inside(m,n):
    return lambda x: np.sin( (pi*m*x)/L ) * (   ( h_bar*h_bar*0.5/m_e )*( n*n*pi*pi/L**2 )*np.sin( n*pi*x/L ) + (a*x/L)*np.sin(n*pi*x/L)  )

# integrates our inside function
def H_mn(m,n):
    return (2/L) * simpson_integral(upper_bound=L, function=H_mn_inside(m,n))

# Analytically gives us values for H with V(x) = ax/L
def H_mn_with_V(m,n):
    if m == n:
        num = h_bar*h_bar * n*n *pi*pi
        den = 2 * m_e * L*L
        return num/den + a/4
    elif m%2 == n%2:
        num = -8*a*m*n
        den = pi*pi * (m*m - n*n)**2
        return num/den
    else:
        return a/2


# =========== Test our Functions ========== #

# We want a 10x10 array of floats for our matrix H
size = 10
H = np.zeros([size, size], float)


# Our Matrix H such that H[S] = E*S
for m in range(size):
    for n in range(size):
        H[m,n] = H_mn_with_V(m+1,n+1)
print(H)

# Calculate our eigenvalues, X
X, V = np.linalg.eigh(H)


X = list( map( lambda x: x/e, X ) )
print("\n\nOur eigenvalues are:\n", X)


