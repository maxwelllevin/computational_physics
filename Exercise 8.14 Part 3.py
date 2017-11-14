# -*- coding: utf-8 -*-
"""
Exercise 8.14 Quantum Oscillators

Created on Tue Nov  7 19:02:33 2017

@author: Maxwell
"""

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
    return v0*x**4/a**4


# Do stuff
def f(r, x, E):
    psi = r[0]
    phi = r[1]
    fpsi = phi
    fphi = (2*m/hbar**2)*(V(x)-E)*psi
    return np.array([fpsi, fphi], float)


# Calculate the wavefunction for a particular energy
def solve(E):
    psi = 0.0
    phi = 1.0
    r = np.array([psi, phi], float)
    for x in np.arange(-L, L, h):
        k1 = h*f(r, x, E)
        k2 = h*f(r+0.5*k1, x+0.5*h, E)
        k3 = h*f(r+0.5*k2, x+0.5*h, E)
        k4 = h*f(r+k3, x+h, E)
        r += (k1+2*k2+2*k3+k4)/6

    return r[0]


def secant(e2):
    # Main program to find the energy using the secant method
    e1 = 0
    psi2 = solve(e1)
    target = e/1000
    while abs(e1-e2) > target:
        psi1, psi2 = psi2, solve(e2)
        e1, e2 = e2, e2-psi2*(e2-e1)/(psi2-psi1)
    # Energy of the Ground State
    print("E =", e2/e, "eV")
    return e2/e


E1 = secant(e)
E2 = secant(400*e)
E3 = secant(910*e)

