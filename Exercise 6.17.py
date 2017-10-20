# -*- coding: utf-8 -*-
"""
Exercise 6.17: Nonlinear Circuits

Created on Thu Oct 19 20:59:32 2017

@author: Maxwell, with help from Iddo Furhmann
"""

import numpy as np


Vp = 5
R1 = 1000
R2 = 4000
R3 = 3000
R4 = 2000
I0 = 3e-9
Vt = 0.05


def f_vec(V1, V2):
    f = np.zeros([2], float)
    f[0] = (V1 - Vp)/R1 + V1/R2 + I0*(np.e**((V1-V2)/Vt) - 1)
    f[1] = (V2 - Vp)/R3 + V2/R4 + I0*(np.e**((V2-V1)/Vt) - 1)
    return f

def jacobian(V1, V2):
    j = np.zeros([2,2], float)
    j[0,0] = 1/R1 + 1/R2 + I0*(np.e**((V1-V2)/Vt))/Vt   # aF0/aV1
    j[0,1] = -I0*(np.e**((V1-V2)/Vt))/Vt                # aF0/aV2
    j[1,0] = -I0*(np.e**((V2-V1)/Vt))/Vt                # aF1/aV1
    j[1,1] = 1/R3 + 1/R4 + I0*(np.e**((V2-V1)/Vt))/Vt   # aF1/aV2
    return j

def mat_inv(mat):
    inv = np.zeros([2, 2], float)
    
    det_inv = 1/(mat[0,0]*mat[1,1] - mat[0,1]*mat[1,0])
    
    inv[0,0] = det_inv*mat[1,1]
    inv[0,1] = -det_inv*mat[0,1]
    inv[1,0] = -det_inv*mat[1,0]
    inv[1,1] = det_inv*mat[0,0]

    return inv





v1 = 0.2
v2 = 1.0

V= np.zeros([2],float)
V[0]=v1
V[1]=v2

for i in range(0,100):
    V -= np.dot( mat_inv(jacobian(V[0],V[1])), f_vec(V[0],V[1]) )
print(V)    