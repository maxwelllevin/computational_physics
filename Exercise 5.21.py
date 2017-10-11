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

# Change x_range and y_range to 50 for part a
x_range = 50
y_range = 50


res = 1
points = 2 * x_range * res
xspace = 2*x_range / (points)
yspace = 2*y_range / (points)

arr = np.zeros([points, points], float)

for i in range(0, points):
    y = i*yspace - y_range
    for j in range(0, points):
        x = j*xspace - x_range
        if (x-5)**2 + y*y == 0 or (x+5)**2 + y*y == 0:
            arr[i,j] = 0
            continue
        arr[i,j] = Phi(x,y)


print("\n\nNegative charge on the left, Positive charge on the right\n")
pl.imshow(arr, origin="lower", extent=[-x_range,x_range, -y_range,y_range])
pl.title("Electric Potential on the X-Y Plane")
pl.xlabel("x, cm")
pl.ylabel("y, cm")
pl.nipy_spectral()
pl.show()


# ======= Make a density plot of E(x,y) (magnitude & direction plots) ====== #

# Magnitude Plot
mag = np.zeros([2*x_range, 2*x_range], float)
for x in range(-x_range, x_range):
    for y in range(-y_range, y_range):
        if (x-5)**2 + y*y == 0 or (x+5)**2 + y*y == 0:
            mag[y,x] = 0
            continue
        field = E_field(x,y)
        mx = field[0]
        my = field[1]
        mag[y+y_range, x+x_range] = np.sqrt(mx*mx + my*my) # Have to reverse the index?

pl.imshow(mag, origin="lower", extent=[-x_range,x_range, -y_range,y_range])
pl.title("Magnitude of E field on the X-Y Plane")
pl.xlabel("x, cm")
pl.ylabel("y, cm")
pl.hsv()
pl.show()


# Direction Plot - HELP
direct = np.zeros([2*x_range, 2*x_range], float)
for x in range(-x_range, x_range):
    for y in range(-y_range, y_range):
        if (x-5)**2 + y*y == 0 or (x+5)**2 + y*y == 0:
            continue
        field = E_field(x,y)
        dirx = field[0]
        diry = field[1]
        direct[y+y_range, x+x_range] = np.arctan(diry/dirx)  #dirx + diry  # Have to reverse the index?

pl.imshow(direct, origin="lower", extent=[-x_range,x_range, -y_range,y_range])
pl.title("Direction of E field on the X-Y Plane")
pl.xlabel("x, cm")
pl.ylabel("y, cm")
pl.show()


# ====== Make a density plot of E field with new charge distribution ======#

print("\nPart c)\n")

#
#
# NOT EVEN CLOSE TO COMPLETE
#
#

L = x_range

# Function that computes the integral of a given function from set bounds with specific depth (~accuracy)
def simpson_integral(lower_bound=0, upper_bound=L, function=lambda x: 1, depth=0.1):
    N = int( (upper_bound - lower_bound) / depth )
    ans = ( function(lower_bound) + function(upper_bound))
    
    for k in range(1, N, 2):
        ans += 4 * function(lower_bound + k * depth)
    
    for k in range(2, N-1, 2):
        ans += 2 * function(lower_bound + k * depth)
    
    ans *= depth / 3
    return ans

def new_E_inside(x1,y1):
    return lambda x,y: (q_0*np.sin(2*pi*x/L)*np.sin(2*pi*y/L) ) / ( 4*pi*e_0*np.sqrt( (x-x1)**2 + (y-y1)**2 ) )

def new_E_x():
    return 0




# Plot of charge distribution
char = np.zeros([2*L, 2*L], float)
for x in range(-L, L):
    for y in range(-L, L):
        char[y+L, x+L] = charge_density(x,y, L, L) # Have to reverse the index?

pl.imshow(char, origin="lower", extent=[-x_range,x_range, -y_range,y_range])
pl.title("Magnitude of E field on the X-Y Plane")
pl.xlabel("x, cm")
pl.ylabel("y, cm")
pl.bone()
pl.show()

