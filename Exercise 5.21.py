# -*- coding: utf-8 -*-
"""
Exrcise 5.21: Electric Field of a Charge Distribution (Pg. 209)

Created on Sat Oct  7 15:16:50 2017

@author: Maxwell, with help from Iddo Fuhrmann 
"""

import numpy as np
import pylab as pl


pi = np.pi
e_0 = 8.85e-12
q = 1
q_0 = 100
L = 5  #L goes from -5 to 5 

# Change x_range and y_range to 50 for part a
x_range = 50
y_range = 50


# =========== FUNCTIONS ============ #


def electric_potential():
    return lambda x,y=0: ( q / ( 4*pi*e_0*np.sqrt( (x-5)**2 + y**2 ) ) ) - ( q / ( 4*pi*e_0*np.sqrt( (x+5)**2 + y**2 ) ) )

def e_x():
    return lambda x,y: -(q/(4*pi*e_0)) * ( (x + 5)/(((x+5)**2 + y*y)**(3/2)) - (x-5)/(((x-5)**2 + y*y)**(3/2))  )

def e_y():
    return lambda x,y: -(q*y/(4*pi*e_0)) * ( ((x+5)**2 + y*y )**(-3/2) - ((x-5)**2 + y*y)**(-3/2) )

def E_field(x,y):
    return ( e_x()(x,y), e_y()(x,y) )

def new_V_inside(x,y,x1,y1):
    return (q_0*np.sin(pi*x/L)*np.sin(pi*y/L) ) / ( 4*pi*e_0*np.sqrt( (x-x1)**2 + (y-y1)**2 ) )

def trapInt(x,y):
    s = 0
    for i in range(-L, L+1):
        for j in range(-L, L+1):
            if(i == L or i == -L):
                c1 = 1/2
            else:
                c1 = 1
            if(j == L or i == -L):
                c2 = 1/2
            else:
                c2 = 1
            if(x != i and y != j):
                s += c1*c2*new_V_inside(i,j, x,y)
    return (s)



# ======= Make a density plot of Phi(x,y) ======= #


Phi = electric_potential()


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
pl.spectral()
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


# Direction Plot
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
#pl.hsv()
pl.show()


# ====== Make a density plot of E field with new charge distribution ======#

print("\nPart c)\n")


points=100
V = np.zeros([points, points], float)
Ex = np.zeros([points, points], float)
Ey = np.zeros([points, points], float)
E = np.zeros([points, points], float)
Edir = np.zeros([points, points], float)

dx=0
dy=0

print(trapInt(0,0))

for i in range(0,100):
    for j in range(0,100):
        Ex[i][j]=(trapInt(i-49.5,j-50)-trapInt(i-50,j-50))/.5
        Ey[i][j]=(trapInt(i-50,j-49.5)-trapInt(i-50,j-50))/.5
        E[i][j]=np.sqrt(Ex[i][j]**2+Ey[i][j]**2)
        Edir[i,j]=np.arctan(Ey[i,j]/Ex[i,j])

pl.imshow(E, origin="lower", extent=[-50,50,-50,50])
pl.bone() 
pl.show()

pl.imshow(Edir, origin="lower", extent=[-50,50,-50,50])
pl.hsv() 
pl.show()