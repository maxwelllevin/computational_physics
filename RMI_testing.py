# -*- coding: utf-8 -*-
"""
RMI test

Created on Sat Sep 30 21:39:20 2017

@author: Maxwell
"""

# HOMEWORK 2 Reading from a data file. Fill in the area for the integrand.
from numpy import loadtxt#, linspace
from pylab import plot, show, xlim, legend, xlabel, ylabel, title
 
 
data = loadtxt('data1.txt')  # loads the data file
T_plot = data[0]  # plot of temperatures
E_1 = data[1]  # plot of first energy
E_2 = data[3]  # plot of second energy
E_3 = data[5]  # plot of third energy
 
count = len(T_plot)  # number of temperature points
 

# Fill in this space with the calculation of RMI for each temperature:
##########      CODE      ##########

print(E_1[0])
print(E_2[0])
print(E_3[0])

T_1 = data[0]
T_2 = data[2]
T_3 = data[4]

print(T_1[1] - T_1[0])
print(T_2[1] - T_2[0])
print(T_3[1] - T_3[0])




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
