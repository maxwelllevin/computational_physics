# -*- coding: utf-8 -*-
"""
Exercise 8.14 Quantum Oscillators

Created on Tue Nov  7 19:02:33 2017

@author: Maxwell
"""

# Repeat all calculations for V(x) = v0x^4 / a^4
# TODO: Calculate normalized wavefunctions for the three states & make a plot



import numpy as np
import pylab as pl


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


def wave(E):
    psi = 0.0
    phi = 1.0
    PSI = []
    r = np.array([psi, phi], float)
    for x in np.arange(-L/2, L/2, h):
        k1 = h*f(r, x, E)
        k2 = h*f(r+0.5*k1, x+0.5*h, E)
        k3 = h*f(r+0.5*k2, x+0.5*h, E)
        k4 = h*f(r+k3, x+h, E)
        r += (k1+2*k2+2*k3+k4)/6
        PSI.append(r[0])
    return np.array(PSI)


def wave_simpson(x_points, y_points):
    N = len(x_points)
    dx = x_points[1]-x_points[0]
    ans = y_points[0] + y_points[-1]
    for k in range(1, N, 2):
        ans += 4 * y_points[k]
    for k in range(2, N-1, 2):
        ans += 2 * y_points[k]
    ans *= dx / 3
    return ans


def secant(e1, e2):
    # Find the energy in a given range using the secant method
    psi2 = solve(e1)
    target = e/1000
    while abs(e1-e2) > target:
        psi1, psi2 = psi2, solve(e2)
        e1, e2 = e2, e2-psi2*(e2-e1)/(psi2-psi1)
    print("E =", e2/e, "eV")
    return e2/e




E0 = secant(0, e)*e
E1 = secant(0, 400*e)*e
E2 = secant(0, 910*e)*e


print(len(wave(E1)))

xPts = np.arange(-L/2, L/2, h)
yPts = wave(E2*e)
norm = wave_simpson(xPts, yPts)
xPts /= norm
yPts /= norm
pl.plot(xPts, yPts)



