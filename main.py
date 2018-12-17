#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 13:37:27 2018

@author: raptev
"""
from math import inf
import optics as o
import numpy as np

import time
"""
    double d0=900e-6;//диаметр накачки, м
    f=1/.254/P*pow(d,2)/pow(d0,2); //qDebug()<<"f"<<f;
    linza(f);
)--l1--FNd--l2--|
"""
start_time = time.clock()

n = 1, 6# nx, ny
l1=20e-3
l2 = 170e-3
rHr = -250e-3
P, y = o.vectorY(*n, 20, 50)
d0 = 900e-6
d = 1e-3

l3 = 1100e-3
f500 = 500e-3
l4 = 1500e-3

HR = o.OptMirror(*n,rHr)
L1 = o.OptSpace( *n, l1)
FNd = o.OptLens(*n, 1/.254/P * d**2/d0**2)
L2 = o.OptSpace(*n, l2)
OM = o.OptMirror(*n, inf, 1)
L3 = o.OptSpace(*n, l3)
F500 = o.OptLens(*n, f500)
L4 = o.OptSpace(*n, l4)

res = HR, L1, FNd, L2, OM
M1 = o.resonatorMatrix(res)
q1 = M1.resonatorMode(532e-9)
opt = res, L3, F500, L4
#print(q1)
#print(q2)

#legend = list(map(lambda y: 'rHr = ' + str(y * 1000) + 'mm', y))
#
#o.grafX(x, q1.rho * 2e3, 'график 1', 'P, W', 'D, mm', legend)

print (time.clock() - start_time, "seconds")
    
#o.graf2d(x,y,q2.rho * 2e6, 'ris1', 'P, Вт', 'радиус зеркала, мм', 150, 200, 1000)
#print (time.clock() - start_time, "seconds")

z, rho, R = o.lineOpt(opt, q1)
legend = list(map(lambda y: 'P = ' + str(y) + 'W', y))
o.grafX(z, rho * 2e3, 'график 1', 'z, mm', 'D, mm',legend)
