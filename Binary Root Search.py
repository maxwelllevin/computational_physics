# -*- coding: utf-8 -*-
"""
Binary Search For Roots

Created on Mon Nov  6 15:56:54 2017

@author: Maxwell Levin
"""

import numpy as np

def f(x):
    return 2 - np.e**(-x)


x1 = -100
x2 = 100
error = 1000
delta = 1e-6

while error > delta:
    xp = (x1 + x2)/2
    if ( f(xp) * f(x1) > 0):
        x1 = xp
    else:
        x2 = xp
    error = np.abs( f(x2) - f(x1))

print((x1 + x2)/2)