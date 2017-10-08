# -*- coding: utf-8 -*-
"""
Exrcise 5.21: Electric Field of a Charge Distribution (Pg. 209)

Created on Sat Oct  7 15:16:50 2017

@author: Maxwell
"""

import numpy as np
import pylab as pl


pi = np.pi
e_0 = 8.85e-12

def electric_potential(q=1):
    return lambda x,y=0: ( q / ( 4*pi*e_0*np.sqrt( (x-5)**2 + y**2 ) ) ) - ( q / ( 4*pi*e_0*np.sqrt( (x+5)**2 + y**2 ) ) )



# ======= Make a density plot of Phi(x,y) ======= #


Phi = electric_potential()

xrange = 15
yrange = 15
points = 2 * xrange
xspace = 2*xrange / (points)
yspace = 2*yrange / (points)

arr = np.empty([points, points], float)

for i in range(0, points):
    y = i*yspace - yrange
    for j in range(0, points):
        x = j*xspace - xrange
        if (x-5)**2 + y*y == 0 or (x+5)**2 + y*y == 0:
            continue
        arr[i,j] = Phi(x,y)



pl.imshow(arr, origin="lower", extent=[-xrange,xrange, -yrange,yrange])
pl.show()


# ======= Make a density plot of E(x,y) (magnitude & direction plots) ====== #






# ======= Fun: Make a plot of Phi at y=0 ======= #
x = np.linspace(-50,50,100)
pl.plot(x, Phi(x))
pl.title("Plot of Electrix Potential Vs x at y=0")
pl.xlabel("x, cm")
pl.ylabel("Phi, V")
pl.show()