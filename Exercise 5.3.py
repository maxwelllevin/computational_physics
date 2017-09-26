# -*- coding: utf-8 -*-
"""
Exercise 5.3: Analytic Integration Example

Created on Fri Sep 22 09:17:36 2017

@author: Maxwell Levin
"""

from numpy import e, linspace
import pylab as pl


def inside_func():
    return lambda t: e**(-t**2)


def simpson_integral(lower_bound, upper_bound, N=100, function=inside_func()):
    depth = (upper_bound - lower_bound) / N
    ans = ( function(lower_bound) + function(upper_bound))
    
    for k in range(1, N, 2):
        ans += 4 * function(lower_bound + k * depth)
    
    for k in range(2, N-1, 2):
        ans += 2 * function(lower_bound + k * depth)
    
    ans *= depth / 3
    return ans




# a) calculate E(x) for 0 <= x <= 3 with a step size of 0.1
print("\na)\n")

step_size = 0.1
our_range = [0,3]
N = int( (our_range[1] - our_range[0]) / step_size )

for i in range(our_range[1] * N + 1):
    f = simpson_integral(our_range[0], our_range[0] + i / N, N)
    f = format(f, '.5f')
    x = format(i / N, '.5f')
    print( "f(", x, ") =", f)


# b) graph E(x) as an actual function:

x = linspace(0, 3, 100)
N = len(x)

pl.plot(x, simpson_integral(0, x, N))
pl.xlabel("X", fontsize=14)
pl.ylabel("E(x)", fontsize=14)
pl.title("Graph of E(x) = integral from 0 to x: e^(-t^2) dt")
pl.show()

