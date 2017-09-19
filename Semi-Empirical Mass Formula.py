# -*- coding: utf-8 -*-
"""
2.10: The Semi-Empirical Mass Formula

Created on Sun Sep 17 14:53:58 2017

@author: Maxwell
"""



# a) take A, Z as input and output binding energy for that atom
print("\nPart a)")

A,Z = int(input("Enter a mass number, A: ")), int(input("Enter an atomic number, Z: "))

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
print("The binding energy is", B, "MeV")





# b) print out the binding energy per nucleon
print("\nPart b)\n")
print("The binding energy per nucleon is:", B / A)





# c) take Z as input and then cycle A through A=Z to A=3Z to find the the largest value of B/A
print("\nPart c)")
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
    

print("The most stable atomic mass number is", A)
print("The binding energy per nucleon would then be", max_BA)
    



    
# d) run Z from 1->100 and print out the most stable value of A for each.
# print out the value of Z for which the binding energy per nucleon is the largest
print("\nPart d)\n")

Z = 1
A = Z
Energy = 0
max_binding_energy = [Z, A, Energy]

# Set a1 -> a4
a1 = 15.8
a2 = 18.3
a3 = 0.714
a4 = 23.2

# Run through Z = 1 to Z= 100 and print most stable A for each
for Z in range(1, 100+1):    
   
    # List to hold our value of A
    local_max = [A, Energy]
   
    for A in range(Z, 3*Z+1):
       
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
        
        # Store the value of Z that yields the most stable value of A overall
        if B / A > max_binding_energy[2]:
            max_binding_energy = [Z, A, B/A]

        # Store the most stable value of A for each Z
        if B / A > local_max[1]:
            local_max = [A, B / A]
            
    print("The most stable A for Z =", Z, "is A =", A)

print("\nThe atomic number that yields the highest binding energy per nucleon is", max_binding_energy[0])
print("Z =", max_binding_energy[0], "yields a binding energy of", max_binding_energy[2], "MeV per nucleon")


