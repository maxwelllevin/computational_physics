# -*- coding: utf-8 -*-
"""
Exercise 2.9: The Madelung constant
TODO: Check if answer is meant to diverge as L->large

Created on Sun Sep 17 11:37:17 2017

@author: Maxwell
"""

e = 1.60217662 * 10**(-19)
perm_vac = 8.85 * 10**(-12)
pi = 3.14159265359
a = 1

def atom_charge(i, j, k):
    """ Returns e if there is a sodium atom, -e if chlorine """
    if i + j + k % 2 == 0: return e
    else: return -e



def atom_potential(i, j, k):
    """ Returns the potential of an atom at i, j, k"""
    numerator = atom_charge(i, j, k)
    denominator = 4 * pi * perm_vac * a * (i**2 + j**2 + k**2)**(0.5)
    return numerator / denominator



def total_potential(L=10):
    """ Returns the sum of all the potentials of atoms in the cubic lattice of length L"""
    
    sum_potential = 0
    
    # Summation
    for i in range(-L, L):
        
        for j in range(-L, L):
            
            for k in range(-L, L):
                
                if i == j == k == 0:
                    # do nothing
                    continue
                sum_potential += atom_potential(i, j, k)
    
    # Return our sum
    return sum_potential


def compute_Madelung(depth=10):
    """ Approximates the Madelung Constant """
    madelung = 4 * pi * perm_vac * a * total_potential(depth) / e
    return madelung



# ====== Test our program ======= #

import datetime
before = datetime.datetime.today()

print (compute_Madelung(30))

after = datetime.datetime.today()
print("Time elapsed = ", after-before)



    