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
q = 1

def electric_potential():
    return lambda x,y=0: ( q / ( 4*pi*e_0*np.sqrt( (x-5)**2 + y**2 ) ) ) - ( q / ( 4*pi*e_0*np.sqrt( (x+5)**2 + y**2 ) ) )

def e_x():
    return lambda x,y: -(q/(4*pi*e_0)) * ( (x + 5)/((x+5)**2 + y*y)**(3/2) - (x-5)/((x-5)**2 + y*y)**(3/2)  )

def e_y():
    return lambda x,y: -(q*y/(4*pi*e_0)) * ( ((x+5)**2 + y*y )**(-3/2) - ((x-5)**2 + y*y)**(-3/2) )

def E_field(x,y):
    return ( e_x()(x,y), e_y()(x,y) )

# ======= Make a density plot of Phi(x,y) ======= #


Phi = electric_potential()

# Change xrange and yrange to 50 for part a
xrange = 20
yrange = 20
points = 2 * xrange
xspace = 2*xrange / (points)
yspace = 2*yrange / (points)

arr = np.zeros([points, points], float)

for i in range(0, points):
    y = i*yspace - yrange
    for j in range(0, points):
        x = j*xspace - xrange
        if (x-5)**2 + y*y == 0 or (x+5)**2 + y*y == 0:
            arr[i,j] = 0
            continue
        arr[i,j] = Phi(x,y)



pl.imshow(arr, origin="lower", extent=[0,2*xrange, 0,2*yrange])
pl.title("Electric Potential on the X-Y Plane")
pl.hsv()
pl.show()


# ======= Make a density plot of E(x,y) (magnitude & direction plots) ====== #

mag = np.zeros([points, points], float)

for i in range(0, points):
    y = i*yspace - yrange
    for j in range(0, points):
        x = j*xspace - xrange
        if (x-5)**2 + y*y == 0 or (x+5)**2 + y*y == 0:
            mag[i,j] = 0
            continue
        field = E_field(x,y)
        mx = field[0]
        my = field[1]
        mag[i,j] = np.sqrt(mx*mx + my*my)

pl.imshow(mag, origin="lower", extent=[0,2*xrange, 0,2*yrange])
pl.title("Magnitude of E field on the X-Y Plane")
pl.bone()
pl.show()






# ======= Fun: Make a plot of Phi at y=0 ======= #
"""
x = np.linspace(-50,50,100)
pl.plot(x, Phi(x))
pl.title("Plot of Electric Potential Vs x at y=0")
pl.xlabel("x, cm")
pl.ylabel("Phi, V")
pl.show()
"""