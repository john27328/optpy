#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 08:07:38 2019

@author: raptev
"""
from math import inf
import optics as o
import numpy as np

def rangeZoom(D,f1,f2,f3):
    f11 = 1/(1/f1+1/f2)
    f12 = f3
    f21 = f1
    f22 = 1/(1/f2 + 1/f3)
    l = f11 + f12, f21 + f22
    if(l[0] < 0 or l[1]<0):
        z=0,0
    else:
        z = D * f12 / f11, D * f22/ f21
    return z

f1 = 55
f2 = -60
f3 = -100
D = 14.4 * 80 / 50

ff1 = np.linspace(40, 60, 100)
ff2 = np.linspace(-100, -30, 100)
ff3 = np.linspace(-100, -20, 100)

z1 = list(map(lambda ff1: rangeZoom(D, ff1, f2,f3),ff1))
z2 = list(map(lambda ff2: rangeZoom(D, f1, ff2,f3),ff2))
z3 = list(map(lambda ff3: rangeZoom(D, f1, f2 ,ff3),ff3))

import matplotlib.pyplot as plt

k = 3
if (k == 1):
    plt.plot(ff1, z1)
elif(k == 2):
    plt.plot(ff2, z2)
elif(k == 3):
    plt.plot(ff3, z3)
        










