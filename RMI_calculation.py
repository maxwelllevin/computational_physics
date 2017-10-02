# -*- coding: utf-8 -*-
"""
RMI_calculation

Created on Fri Sep 22 10:04:10 2017

@author: Maxwell Levin 
"""
# HOMEWORK 2 Reading from a data file. Fill in the area for the integrand.
from numpy import loadtxt
from pylab import plot, show, xlim, legend, xlabel, ylabel, title
 
 
data = loadtxt('data1.txt')  # loads the data file
T_plot = data[0]  # plot of temperatures
E_1 = data[1]  # plot of first energy
E_2 = data[3]  # plot of second energy
E_3 = data[5]  # plot of third energy
 
count = len(T_plot)  # number of temperature points
 

# Fill in this space with the calculation of RMI for each temperature:
##########      CODE      ##########

# Contains temperature data
T_data = data[0]

# Get Tmin and Tmax from the data
T_min = T_data[0]
T_max = T_data[-1]

# Find the length of our temperature data array
T_length = len( T_data )


# The inside function we'll use as our integrand
def inside_RMI():
    return lambda i: (2*E_1[i] - E_2[i] - 2*E_3[i]) / (T_data[i]**2)

# Computes the integral of the integrand using Simpson's Rule
def simpson_RMI(index_s, index_f, func_length, func=lambda x:1):
    integral = func(index_s) + func(index_f)
    for k in range(1, func_length, 2):
        integral += 4 * func(index_s + k)
    for k in range(2, func_length, 2):
        integral += 2 * func(index_s + k)
    return integral / 3


# Builds our RMI function of 1-variable (integer input only)
def RMI_builder():
    return lambda index: simpson_RMI(index, -1, T_length - index, inside_RMI())

# Our RMI function (integer input from 0 to 997 only)
RMI = RMI_builder()

# Test RMI at a specific index
#print( RMI(0) )

# Calculate RMI at all the indices
RMI_plot = []
for i in range(T_length):
    RMI_plot.append( RMI(i) )

# Plot our results
plot(T_plot, RMI_plot, label='RMI')
xlabel('Temperature (K)')
ylabel('RMI (J)')
title('RMI')
xlim(0,100)
show()

# Error for our Simpson's rule is h^4, our h = 0.1
print("Our estimated error is", format(0.1**4, '0.4f'), "at each point we've calculated RMI(T)")


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