# -*- coding: utf-8 -*-
"""
Adaptive Method for Integration with Trapezoidal Sums

Created on Wed Sep 27 15:19:45 2017

@author: Maxwell Levin
"""

import numpy as np


def trapezoidal_error(function=lambda x:1, delta=1e-5, lower=0, upper=1, h=0.1):
    
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
    
integral, error = trapezoidal_error(lambda x:x**2, 1e-9)
print(integral, error)