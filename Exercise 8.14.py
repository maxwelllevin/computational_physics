# -*- coding: utf-8 -*-
"""
Exercise 8.14 Asymmetric Quantum Well

Created on Tue Nov  7 19:02:33 2017

@author: Maxwell
"""

# TODO: Calculate Energies of first two excited states
# TODO: Repeat all calculations for new V(x) = v0x^4 / a^4
# TODO: Calculate normalized wavefunctions for the three states & make a plot

import numpy as np

# Constants
m = 9.1094e-31
hbar = 1.0546e-34
e = 1.6022e-19
L = 1e-10
N = 1000
h = L/N
v0 = 50*e
a = 1e-11

# Potential Function
def V(x):
    return v0*x*x/a**2

def f(r,x,E):
    psi = r[0]
    phi = r[1]
    fpsi = phi
    fphi = (2*m/hbar**2)*(V(x)-E)*psi
    return np.array([fpsi, fphi], float)

# Calculate the wavefunction for a particular energy
def solve(E):
    psi = 0.0
    phi = 1.0
    r = np.array([psi,phi],float)
    
    for x in np.arange(-L,L,h):
        k1 = h*f(r,x,E)
        k2 = h*f(r+0.5*k1, x+0.5*h, E)
        k3 = h*f(r+0.5*k2, x+0.5*h, E)
        k4 = h*f(r+k3,x+h,E)
        r += (k1+2*k2+2*k3+k4)/6
    
    return r[0]

# Main program to find the energy using the secant method
E1 = 0.0
E2 = e
psi2 = solve(E1)

target = e/1000
while abs(E1-E2) > target:
    psi1,psi2 = psi2, solve(E2)
    E1,E2 = E2,E2-psi2*(E2-E1)/(psi2-psi1)

# Energy of the Ground State
print("E =", E2/e, "eV")