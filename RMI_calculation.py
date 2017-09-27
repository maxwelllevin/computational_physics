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
 

Tmax = 997
Tmin = 0
h = 0.1
 
print(E_1[-1])
 
# Fill in this space with the calculation of RMI for each temperature:
##########      CODE      ##########
 

# Returns the value of E1 at an index
def get_E1(T):
    index = int( (Tmax - T) / h )
    return E_1[index]

# Returns the value of E2 at an index
def get_E2(T):
    index = int( (Tmax - T) / h )
    return E_2[index]

# Returns the value of E3 at an index
def get_E3(T):
    index = int( (Tmax - T) / h )
    return E_3[index]


# Builds our inside function of T
def inside_func():
    return lambda T: (2*get_E1(int(T)) - get_E2(int(T)) - 2 * get_E3(int(T))) / (T**2)


# Our integral computation function:
def simpson_integral(lower_bound, upper_bound, N=100, function=inside_func()):
    depth = (upper_bound - lower_bound) / N
    ans = ( function(lower_bound) + function(upper_bound))
    
    for k in range(1, N, 2):
        ans += 4 * function(lower_bound + k * depth)
    
    for k in range(2, N-1, 2):
        ans += 2 * function(lower_bound + k * depth)
    
    ans *= depth / 3
    return ans






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