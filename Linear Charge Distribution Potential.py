# -*- coding: utf-8 -*-
"""
HW#2: Potential & Electric Field of Linear Charge Distribution

Created on Tue Sep 26 14:23:01 2017

@author: Maxwell Levin
"""

import numpy as np
import pylab as pl


# Inside function for the potential integral
def V_inside(x, y, z):
    return lambda theta: (1 / (4*np.pi)) * np.sqrt( 4*np.cos(theta)**2 + 9*np.sin(theta)**2 ) / np.sqrt( (x - 3*np.cos(theta))**2 + (y - 2*np.sin(theta))**2 + z**2 )

# Inside function for x-component of E field
def E_x_inside(x, y, z):
    return lambda theta: (1 / (4*np.pi)) * ( x - 3*np.cos(theta) ) * np.sqrt( 4*np.cos(theta)**2 + 9*np.sin(theta)**2 ) / ( (x - 3*np.cos(theta))**2 + (y - 2*np.sin(theta))**2 + z**2 )**1.5

# Inside function for y-component of E field
def E_y_inside(x, y, z):
    return lambda theta: (1 / (4*np.pi)) * ( y - 2*np.sin(theta) ) * np.sqrt( 4*np.cos(theta)**2 + 9*np.sin(theta)**2 ) / ( (x - 3*np.cos(theta))**2 + (y - 2*np.sin(theta))**2 + z**2 )**1.5

# Inside function for z-component of E field
def E_z_inside(x, y, z):
    return lambda theta: (1 / (4*np.pi)) * z * np.sqrt( 4*np.cos(theta)**2 + 9*np.sin(theta)**2 ) / ( (x - 3*np.cos(theta))**2 + (y - 2*np.sin(theta))**2 + z**2 )**1.5


# Function that computes the integral of a given function from set bounds with specific depth (~accuracy)
def simpson_integral(lower_bound, upper_bound, function=lambda x: 1, depth=0.1):
    #depth = (upper_bound - lower_bound) / N
    N = int( (upper_bound - lower_bound) / depth )
    ans = ( function(lower_bound) + function(upper_bound))
    
    for k in range(1, N, 2):
        ans += 4 * function(lower_bound + k * depth)
    
    for k in range(2, N-1, 2):
        ans += 2 * function(lower_bound + k * depth)
    
    ans *= depth / 3
    return ans

# Function that computes the integral of a given function from set bounds with specific starting depth and accuracy requirement of result
def trapezoidal_error(lower=0, upper=1, function=lambda x: 1, h=0.1, delta=1e-7):
    
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
    
    # Return our integral and the last computed error
    return I_1, error


# ============ Find the V and E at (x,y,z), and print the error of each =========== #
print("\nPart a)")
x,y,z = 1,4,7
#x,y,z = int(input( "Enter an x,y,z: ")), int(input()), int(input())

# Create functions for our integrals
V = V_inside(x,y,z)
Ex = E_x_inside(x,y,z)
Ey = E_y_inside(x,y,z)
Ez = E_z_inside(x,y,z)

# set our depth to be 0.1
h = 0.1

# Compute the Potential function at x,y,z
ansV = simpson_integral(0, 2*np.pi, V, h)

# Compute the error 
errorV = h**4

# Print our results
print("V at", (x,y,z), "is:", ansV)
print("Error of V is", format(errorV, '0.6f') )

# Calculate the E field at x,y,z
ansEx = simpson_integral(0, 2*np.pi, Ex, h)
ansEy = simpson_integral(0, 2*np.pi, Ey, h)
ansEz = simpson_integral(0, 2*np.pi, Ez, h)

# Calculate the error of Simpson's rule (h^4 always)
errorE = format(h**4, '.6f')

# Print our results
print("E at", (x,y,z), "is: <", ansEx,",", ansEy,",", ansEz,">")
print("The error of the E field is: <", errorE, ",", errorE,",", errorE, ">")




# ========== Use the adaptive trapezoidal method to obtain 6-digit accuracy for V and E ========== #
print("\nPart b)")
x,y,z = 1,2,5
#x,y,z = int(input( "Enter an x,y,z: ")), int(input()), int(input())

# Create functions for our integrals
V = V_inside(x,y,z)
Ex = E_x_inside(x,y,z)
Ey = E_y_inside(x,y,z)
Ez = E_z_inside(x,y,z)

# set our depth to be 0.1
h = 0.1

# set our accuracy to be 1e-6
accuracy = 1e-6

# Calculate and print our Potential along with its error
ansV = trapezoidal_error(0, 2*np.pi, V, h, accuracy)
print("V at", (x,y,z), "is:", ansV[0], "with error:", ansV[1])

# Calculate the E field at x,y,z
ansEx = trapezoidal_error(0, 2*np.pi, Ex, h, accuracy) 
ansEy = trapezoidal_error(0, 2*np.pi, Ey, h, accuracy)
ansEz = trapezoidal_error(0, 2*np.pi, Ez, h, accuracy)

# Print our E field along with its error
print("E at", (x,y,z), "is: <", ansEx[0], ",", ansEy[0],",", ansEz[0], ">")
print("The error of E is: <", ansEx[1], ",", ansEy[1], ",", ansEz[1], ">")




# ======== Make a density plot ========= #
print("\nPart c)")

# Make our potential function for our density plot
def V_potential(x,y):
    return simpson_integral(0, 2*np.pi, V_inside(x,y,1))

xrange = 10
yrange = 10
points = 100
xspace = xrange / (2*points)
yspace = yrange / (2*points)

arr = np.empty([points, points], float)

for i in range(-points, points):
    y = yspace*i
    for j in range(-points, points):
        x = xspace*j
        arr[i,j] = V_potential(x,y)

pl.imshow(arr, origin="lower", extent=[-xrange,xrange,-yrange,yrange])
pl.bone() 
pl.show()
















