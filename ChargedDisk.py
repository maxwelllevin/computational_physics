# -*- coding: utf-8 -*-
"""
Physics 251 Electric Field of Charged Disk

Created on Wed Sep 20 20:43:47 2017

@author: Maxwell
"""

import numpy as np
import pylab as pl


sigma = 0.5 * 10**(-9)
eta0 = 8.85 * 10**(-12)

def build_field_function(radius):
    return lambda z: (sigma/(2*eta0))*(1 - (1 + radius**2 / z**2)**(-0.5))

num_19a = build_field_function(.3)

x = np.linspace(0.01, 1, 100)

pl.plot(x, num_19a(x), '-b', label='E_exact')
pl.plot(x, x*sigma/(2*eta0*x), '-r', label='E_approx')


pl.title("Electric Field Strength, E, versus Distance, X")
pl.xlabel("X (m)")
pl.ylabel("E (C / N)")
pl.legend()
pl.savefig("num19a.pdf")