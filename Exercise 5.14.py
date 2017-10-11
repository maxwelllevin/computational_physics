# -*- coding: utf-8 -*-
"""
Exercise 5.14

Created on Sun Oct  8 19:25:17 2017

@author: Maxwell, with help from Evan Davidson and Eli Barnes
"""

import numpy as np
import pylab as pl

# Declare our constants
# a --> d are our bounds of integration / size of the plate
a = -5
b = 5
c = -5
d = 5
G = 6.674e-11
sigma = 100

# Number of iterations
N = 100

# Function inside our double integral
def func():
    return lambda x,y,z: (x**2+y**2+z**2)**(-3/2)

# Function to compute 2-dimensional integration with the Trapezoidal Rule
def trapezoidal_2D(z, n=N, f = func()) :
    h = (b - a) / n
    k = (d - c) / n
    I = f(a, c, z) + f(a, c, z) + f(a, d, z) + f(b, d, z)
    
    for i in range(n) :
        xi = a+ i*h
        yj = c +i*k
        I += 2 * ( f(xi, c, z) + f(xi, d, z) + f(a, yj, z) + f(b, yj, z) )
    
    for i in range(n):
        for j in range(n) :
            I += 4 * f(a+i*h, c+j*k, z)
    
    # Multiply our integral by constants and return it        
    I *= G * z * h * k * sigma / 4
    return I


# ======= Make our plot ======= #

y = []
x = np.arange(.1,10.1,.1)



for i in range(1,101) :
    z = i*.1
    fz = trapezoidal_2D(z)
    y.append(fz)


pl.plot(x,y)
pl.xlabel("Z, m")
pl.ylabel("Fz, N")
pl.title("Plot of Force on a Point Mass")
pl.show()