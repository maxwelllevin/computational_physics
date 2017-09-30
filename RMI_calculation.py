# -*- coding: utf-8 -*-
"""
RMI_calculation

Created on Fri Sep 22 10:04:10 2017

@author: Maxwell Levin
"""

# HOMEWORK 2 Reading from a data file. Fill in the area for the integrand.
from numpy import loadtxt, linspace
from pylab import plot, show, xlim, ylim, legend, xlabel, ylabel, title
 
 
data = loadtxt('data1.txt')  # loads the data file
T_plot = data[0]  # plot of temperatures
E_1 = data[1]  # plot of first energy
E_2 = data[3]  # plot of second energy
E_3 = data[5]  # plot of third energy
 
count = len(T_plot)  # number of temperature points
 

# Fill in this space with the calculation of RMI for each temperature:
##########      CODE      ##########
 

"""

def simpson_integral(lower_bound, upper_bound, N=100, function=inside_func()):
    depth = (upper_bound - lower_bound) / N
    ans = ( function(lower_bound) + function(upper_bound))
    
    for k in range(1, N, 2):
        ans += 4 * function(lower_bound + k * depth)
    
    for k in range(2, N-1, 2):
        ans += 2 * function(lower_bound + k * depth)
    
    ans *= depth / 3
    return ans

"""

# Make sure data is printing correctly
#print( E_1[0] )
#print( E_2[0] )
#print( E_3[0] )

# Initialize Constants
Tmin = 0.01
Tmax = 100
deltaT = 0.1

# Initialize Variables : SOMETHING BROKEN (DEBUG)
T = 100
index = int( 10 * (T - Tmin) )
print("Index starts at:", index)
integral = ( 2*E_1[0] - E_2[0] - 2*E_3[0] ) / Tmin**2  + ( 2*E_1[-1] - E_2[-1] - 2*E_3[-1] ) / Tmax**2

# Compute RMI(T) for any T
while T < Tmax and 0 <= index < len(E_1):
    
    if index % 2 == 0:
        integral += 4*( 2*E_1[index] - E_2[index] - 2*E_3[index] ) / T**2
    else:
        integral += 2*( 2*E_1[index] - E_2[index] - 2*E_3[index] ) / T**2
    
    integral *= deltaT/3
    T += deltaT
    index += 1

# PRINTS -1.05192287916e-06 IF T<100 AND 10093.0139622 IF T>=100.. WHY???
print("Index got to:", index)
print("T got to:", T)
print("The integral is:", integral)



def x_squared(x):
    return x**2
arr = [1, 2, 3, 4, 5]

mapped = map(x_squared, arr)
mapped = list(mapped)

print(mapped)


##########      CODE      ##########
 
plot(T_plot, E_1, label='Energy 1')
plot(T_plot, E_2, label='Energy 2')
plot(T_plot, E_3, label='Energy 3')
xlabel('Temperature (K)')
ylabel('Energy (J)')
title('Energy')
xlim(0, 100)
legend()
show()