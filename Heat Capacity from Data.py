# -*- coding: utf-8 -*-
"""
Heat Capacity from Data

Created on Wed Oct 11 20:17:35 2017

@author: Maxwell
"""

from numpy import loadtxt
from pylab import plot, show, xlim, legend, xlabel, ylabel, title
 
 
data = loadtxt('data1.txt')  # loads the data file
T_plot = data[0]  # plot of temperatures
E_1 = data[1]  # plot of first energy
E_2 = data[3]  # plot of second energy
E_3 = data[5]  # plot of third energy
 
count = len(T_plot)  # number of temperature points
 

##########      CODE      ##########

# Find the length of our temperature data array
T_length = len( T_plot )

# List to store our heat capacity
# C = d/dT[E_3(T)]
HC = []

# Our step size for Temperature (0.1)
h = T_plot[1] - T_plot[0]

# Calculate the derivative
for i in range(T_length - 1):
    # Forward Difference
    tempC = E_3[i+1] - E_3[i]
    tempC /= h
    HC.append(tempC)
# Backward Difference for the last element
tempC = E_3[-1] - E_3[-2]
tempC /= h
HC.append(tempC)


# Plot our Heat Capacity Data
plot(T_plot, HC, 'r')
xlabel('Temperature (K)')
ylabel('Heat Capacity (J)')
title('Heat Capacity of E_3')
xlim(0, 100)
show()
    


##########      CODE      ##########
 
#plot(T_plot, E_1, label='Energy 1')
#plot(T_plot, E_2, label='Energy 2')
plot(T_plot, E_3, label='Energy 3')
xlabel('Temperature (K)')
ylabel('Energy (J)')
title('Energy')
xlim(0, 100)
legend()
show()