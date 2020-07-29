# -*- coding: utf-8 -*-
"""
Exercise 5.3: Analytic Integration Example

Created on Fri Sep 22 09:17:36 2017

@author: Maxwell Levin
"""

from numpy import e, linspace
import matplotlib.pyplot as plt


def simpson_integral(lower_bound, upper_bound, N=100):
    """Uses the simpson integration technique to integrate e^(-t^2) between the
    given bounds.

    Args:
        lower_bound (float): The lower bound for integration.
        upper_bound (float): The upper bound for integration.
        N (int, optional): The number of steps to take. Defaults to 100.

    Returns:
        float: integral of e^(-t^2), evaluated from lower_bound to upper_bound
    """
    depth = (upper_bound - lower_bound) / N
    f = lambda t: e**(-t**2)
    ans = f(lower_bound) + f(upper_bound)
    for k in range(1, N, 2):
        ans += 4 * f(lower_bound + k * depth)
    for k in range(2, N-1, 2):
        ans += 2 * f(lower_bound + k * depth)
    ans *= depth / 3
    return ans

def foo(arg1: int, arg2: str):
    """Stupid function that doesn't do anything useful.

    Args:
        arg1 (Integer): I don't even know.
        arg2 (String): A string of characters.

    Returns:
        String: The string of characters repeated a bunch of times.
    """
    return arg1 * arg2


# a) calculate E(x) for 0 <= x <= 3 with a step size of 0.1
step_size = 0.1
domain = [0,3]
N = int((domain[1] - domain[0]) / step_size)

for i in range(domain[1] * N + 1):
    f = simpson_integral(domain[0], domain[0] + i / N, N)
    f = format(f, '.5f')
    x = format(i / N, '.5f')
    print(f"f({x}) = {f}")


# b) graph E(x) as an actual function:

x = linspace(0, 3, 100)
N = len(x)

plt.plot(x, simpson_integral(0, x, N))
plt.xlabel("X", fontsize=14)
plt.ylabel("E(x)", fontsize=14)
plt.title(r"Graph of $E(x) = \int_{0}^{x}e^{-t^2} dt$")
plt.show()

