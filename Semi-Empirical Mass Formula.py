# -*- coding: utf-8 -*-
"""
2.10: The Semi-Empirical Mass Formula

Created on Sun Sep 17 14:53:58 2017

@author: Maxwell
"""



# a) take A, Z as input and output binding energy for that atom
"""
A,Z = int(input("Enter a value for A: ")), int(input("Enter a value for Z: "))

# Set a1 -> a4
a1 = 15.8
a2 = 18.3
a3 = 0.714
a4 = 23.2

# Set a5
if A % 2 == 0:
    if Z % 2 == 0:
        a5 = 12.0
    else:
        a5 = -12
else:
    a5 = 0


# Calculate binding energy
B = a1 * A - a2 * A**(2/3) - a3 * (Z**2) / (A**(1/3)) - a4 * ((A - 2 * Z)**2) / A + a5 / (A**(1/2))

# print the result
print("The binding energy is:", B, "MeV")
"""




"""
# b) print out the binding energy per nucleon

print("The binding energy per nucleon is:", B / A)
"""




# c) take Z as input and then cycle A through A=Z to A=3Z to find the the largest value of B/A
#

Z = int(input("Enter an atomic number, Z: "))

# Set a1 -> a4
a1 = 15.8
a2 = 18.3
a3 = 0.714
a4 = 23.2

# Cycle through A
max_BA = 0
for A in range(Z, 3 * Z + 1):
    # Set a5
    if A % 2 == 0:
        if Z % 2 == 0:
            a5 = 12.0
        else:
            a5 = -12
    else:
        a5 = 0
    
    # Calculate binding energy
    B = a1 * A - a2 * A**(2/3) - a3 * (Z**2) / (A**(1/3)) - a4 * ((A - 2 * Z)**2) / A + a5 / (A**(1/2))
    
    temp_BA = B / A
    if temp_BA > max_BA:
        max_BA = temp_BA
        A_in_BA = A
    

print("A =", A)
print("Binding energy per nucleon =", max_BA)
    

    
    
