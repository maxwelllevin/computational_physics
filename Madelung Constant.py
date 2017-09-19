# -*- coding: utf-8 -*-
"""
Exercise 2.9: The Madelung constant

Created on Sun Sep 17 11:37:17 2017

@author: Maxwell
"""

# NOTE: This program is slow, but does not use much memory. 
# May be able to increase speed by predetermining -1 ^ i+j+k,
# manually accounting for i=j=k=0, and analyzing symmetries



# get L input
L = int( input("Enter a value for L: ") )

# initialize madelung
M = 0


# Loop through all integer lattice points in a cube of length 'L'
# compute the contribution of each atom at each point and add to 'M'
i = -L
while i < L+1:
    j = -L
    while j < L+1:
        k = -L
        while k < L+1:
            # Ignore the origin
            if (i==0) and (j==0) and (k==0):
                k += 1
            else:
                # If i+j+k even, atom is Na+, else Cl-
                # Divide by the distance from current point to NaCl at origin
                M += ((-1)**(i+j+k)) / ( i**2 + j**2 + k**2 )**0.5
                k += 1
        j += 1
    i += 1


# Print the result
print("The MadeLung Constant is approximately", M)
