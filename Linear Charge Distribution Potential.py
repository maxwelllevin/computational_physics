# -*- coding: utf-8 -*-
"""
HW#2: Potential & Electric Field of Linear Charge Distribution

Created on Tue Sep 26 14:23:01 2017

@author: Maxwell Levin
"""

import numpy as np
import pylab as pl

def V_inside(x, y, z):
    return lambda theta: (1 / (4*np.pi)) * np.sqrt( 4*np.cos(theta)**2 + 9*np.sin(theta)**2 ) / np.sqrt( (x - 3*np.cos(theta))**2 + (y - 2*np.sin(theta))**2 + z**2 )

def E_inside(x, y, z):
    return lambda theta: (1 / (4*np.pi)) * np.sqrt( 4*np.cos(theta)**2 + 9*np.sin(theta)**2 ) / ( (x - 3*np.cos(theta))**2 + (y - 2*np.sin(theta))**2 + z**2 )**1.5


def simpson_integral(lower_bound, upper_bound, depth=0.1, function=lambda x: 1):
    #depth = (upper_bound - lower_bound) / N
    N = int( (upper_bound - lower_bound) / depth )
    ans = ( function(lower_bound) + function(upper_bound))
    
    for k in range(1, N, 2):
        ans += 4 * function(lower_bound + k * depth)
    
    for k in range(2, N-1, 2):
        ans += 2 * function(lower_bound + k * depth)
    
    ans *= depth / 3
    return ans



# ============ Test: find the value of V at (x,y,z) = user input(x,y,z) =========== #
x,y,z = 1,4,7
#x,y,z = int(input( "Enter an x,y,z: ")), int(input()), int(input())
V = V_inside(x,y,z)
E = E_inside(x,y,z)

print(V(2*np.pi))


# set our depth to be 0.1
h = 0.1

# Our answer is the simpson integral of our function from 0 to 2pi
ansV = simpson_integral(0, 2*np.pi, h, V)
print("V at", (x,y,z), "is:", ansV)

ansE = simpson_integral(0, 2*np.pi, h, E)
print("E at", (x,y,z), "is:", ansE)



# ========== TODO: use adaptive trapezoidal method to obtain 6-digit accuracy ========== #







# ======== TODO: make a density plot ========= #



