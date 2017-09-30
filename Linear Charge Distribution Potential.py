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

# TODO: Fix the numerator (r - r') to be vector subraction, not distance
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


def trapezoidal_error(lower=0, upper=1, h=0.1, delta=1e-4, function=lambda x: 1):
    
    # Set N to be the number of steps needed
    N = int( (upper - lower) / h )
    
    # Calculate I_1
    I_1 = ( 0.5*function(lower) + 0.5*function(upper) )
    for k in range(1,N):
        I_1 += function(lower + k*h )
    I_1 *= h
    
    # Declare our error tracker
    error = 10
    
    # Loop while our error is undesirable
    while error > delta:
        h /= 2
        N *= 2
        
        S = 0
        for k in range(1, N, 2):
            S += function(lower + k*h )
        I_2 = 0.5*I_1 + h*S
        
        
        error = (1/3) * (I_2 - I_1)
        error = np.abs(error)
        
        I_1 = I_2
    
    # Return our integral 
    return I_1, error


# ============ Test: find the value of V at (x,y,z) = user input(x,y,z) =========== #
print("\nPart a)")
x,y,z = 1,4,7
#x,y,z = int(input( "Enter an x,y,z: ")), int(input()), int(input())
V = V_inside(x,y,z)
E = E_inside(x,y,z)

# set our depth to be 0.1
h = 0.1

# Our answer is the simpson integral of our function from 0 to 2pi
ansV = simpson_integral(0, 2*np.pi, h, V)
print("V at", (x,y,z), "is:", ansV)

ansE = simpson_integral(0, 2*np.pi, h, E)
print("E at", (x,y,z), "is:", ansE)



# ========== TODO: use adaptive trapezoidal method to obtain 6-digit accuracy ========== #
print("\nPart b)")
x,y,z = 1,4,7
#x,y,z = int(input( "Enter an x,y,z: ")), int(input()), int(input())
V = V_inside(x,y,z)
E = E_inside(x,y,z)

# set our depth to be 0.1
h = 0.1

# set our accuracy to be 1e-6
accuracy = 1e-6

# Our answer is the simpson integral of our function from 0 to 2pi
ansV = trapezoidal_error(0, 2*np.pi, h, accuracy, V)
print("V at", (x,y,z), "is:", ansV[0], "with error:", ansV[1])

ansE = trapezoidal_error(0, 2*np.pi, h, accuracy, E)
print("E at", (x,y,z), "is:", ansE[0], "with error:", ansE[1])






# ======== TODO: make a density plot ========= #



