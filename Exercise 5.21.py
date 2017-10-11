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
q_0 = 100 

def electric_potential():
    return lambda x,y=0: ( q / ( 4*pi*e_0*np.sqrt( (x-5)**2 + y**2 ) ) ) - ( q / ( 4*pi*e_0*np.sqrt( (x+5)**2 + y**2 ) ) )

def e_x():
    return lambda x,y: -(q/(4*pi*e_0)) * ( (x + 5)/(((x+5)**2 + y*y)**(3/2)) - (x-5)/(((x-5)**2 + y*y)**(3/2))  )
    #return lambda x,y: -(q/(4*pi*e_0))*((x+5)/(((x+5)**(2)+y**2)**(3/2)))-(x-5)/(((x-5)**(2)+y**2)**(3/2))

def e_y():
    return lambda x,y: -(q*y/(4*pi*e_0)) * ( ((x+5)**2 + y*y )**(-3/2) - ((x-5)**2 + y*y)**(-3/2) )
    #return lambda x,y: -(q/(4*pi*e_0))*(y/(((x+5)**(2)+y**2)**(3/2))-y/(((x-5)**(2)+y**2)**(3/2)))

def E_field(x,y):
    return ( e_x()(x,y), e_y()(x,y) )

def charge_density(x,y, Lx=5, Ly=5):
    return q_0 * np.sin( 2*pi*x / Lx ) * np.sin( 2*pi*y / Ly ) 


# ======= Make a density plot of Phi(x,y) ======= #


Phi = electric_potential()

# Change xrange and yrange to 50 for part a
xrange = 15
yrange = 15


res = 1
points = 2 * xrange * res
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


print("\n\nNegative charge on the left, Positive charge on the right\n")
pl.imshow(arr, origin="lower", extent=[-xrange,xrange, -yrange,yrange])
pl.title("Electric Potential on the X-Y Plane")
pl.xlabel("x, cm")
pl.ylabel("y, cm")
pl.nipy_spectral()
pl.show()


# ======= Make a density plot of E(x,y) (magnitude & direction plots) ====== #

# Magnitude Plot
mag = np.zeros([2*xrange, 2*xrange], float)
for x in range(-xrange, xrange):
    for y in range(-yrange, yrange):
        if (x-5)**2 + y*y == 0 or (x+5)**2 + y*y == 0:
            mag[y,x] = 0
            continue
        field = E_field(x,y)
        mx = field[0]
        my = field[1]
        mag[y+yrange, x+xrange] = np.sqrt(mx*mx + my*my) # Have to reverse the index?

pl.imshow(mag, origin="lower", extent=[-xrange,xrange, -yrange,yrange])
pl.title("Magnitude of E field on the X-Y Plane")
pl.xlabel("x, cm")
pl.ylabel("y, cm")
pl.hsv()
pl.show()


# Direction Plot - HELP
direct = np.zeros([2*xrange, 2*xrange], float)
for x in range(-xrange, xrange):
    for y in range(-yrange, yrange):
        if (x-5)**2 + y*y == 0 or (x+5)**2 + y*y == 0:
            direct[y,x] = 0
            continue
        field = E_field(x,y)
        dirx = field[0]
        diry = field[1]
        direct[y+yrange, x+xrange] = dirx + diry  # Have to reverse the index?

pl.imshow(direct, origin="lower", extent=[-xrange,xrange, -yrange,yrange])
pl.title("Direction of E field on the X-Y Plane")
pl.xlabel("x, cm")
pl.ylabel("y, cm")
pl.show()


# ====== Make a density plot of E field with new charge distribution ======#


#
#
# NOT EVEN CLOSE TO COMPLETE
#
#

Lx = xrange
Ly = yrange

# Plot of charge distribution
char = np.zeros([2*Lx, 2*Ly], float)
for x in range(-Lx, Lx):
    for y in range(-Ly, Ly):
        char[y+Ly, x+Lx] = charge_density(x,y, Lx, Ly) # Have to reverse the index?

pl.imshow(char, origin="lower", extent=[-xrange,xrange, -yrange,yrange])
pl.title("Magnitude of E field on the X-Y Plane")
pl.xlabel("x, cm")
pl.ylabel("y, cm")
pl.bone()
pl.show()

