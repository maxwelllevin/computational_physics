# -*- coding: utf-8 -*-
"""
Exercise 6.16: The Lagrange Point
Created on Wed Oct 18 21:20:29 2017
@author: Maxwell
"""

""" 
We seek the root of GM/r^2 -Gm/(R-r)^2 = (w^2)r which lies in between 
the position of the Earth and the Moon. We use Newton's method to 
calculate the Lagrange point (the root) to four decimal places.
"""


import numpy as np
import pylab as pl


G = 6.674e-11
M = 5.974e+24
m = 7.348e+22
R = 3.844e+8
w = 2.662e-6


def F(r):
    return G*M*((R-r)**2) - G*m*r*r - w*w*((R-r)**2)*(r**3)

def LHS(r):
    return G*M*((R-r)**2)

def RHS(r):
    return G*m*r*r + w*w*((R-r)**2)*(r**3)
    


# This is just for visualization purposes
domain = np.arange(0,R,1e+6)
f = []
g = []
for i in domain:
    f.append(LHS(i))
    g.append(RHS(i))
pl.plot(domain,f)
pl.plot(domain,g)
pl.show()




r1 = 3.4e+8     # my guess based on the graph
h = 1e-4        # used for calculating the derivative

for i in range(1000):
    fr = F(r1)                  # value of the function at r1
    df = (F(r1+h) - F(r1-h))/h  # value of the derivative at r1
    r1 -= fr/df                 # update the value of r1

r1 /= 1e+8
print("The Lagrange Point is at: %.5f" %r1, "e+8 m")            # print r1