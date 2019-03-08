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

n = 1, 20# nx, ny
f1 = 40E-2
f2 = 20E-2
f3 = -50e-2
y1 = o.OptY(*n, 14.3e-2 / 2, 0)
l1 = 100e-2
L1 = o.OptSpace(*n,l1);

F1 = o.OptLens(*n, f1)
l2, y = o.vectorY(*n, 81, 110)
L2 = o.OptSpace(*n, l2 * 1e-2)
F2 = o.OptLens(*n,f2)

opt = L1, F1, L2, F2
M = o.linMatrix(opt)
l3 = f3 + (M * y1).focus()
L3 = o.OptSpace(*n, l3)
F3 = o.OptLens(*n, f3)
opt2 = opt, L3, F3, L1
#print(q1)
#print(q2)

#legend = list(map(lambda y: 'rHr = ' + str(y * 1000) + 'mm', y))
#
#o.grafX(x, q1.rho * 2e3, 'график 1', 'P, W', 'D, mm', legend)

print (time.clock() - start_time, "seconds")
    
#o.graf2d(x,y,q2.rho * 2e6, 'ris1', 'P, Вт', 'радиус зеркала, мм', 150, 200, 1000)
#print (time.clock() - start_time, "seconds")

z, rho = o.lineOpt(opt2, y1)
if (1):
    legend = list(map(lambda y,l: 'L1 = ' + str(y) + 'mm' + ' L2 = ' + str(l * 100) + 'mm', y,l3[0]))
    o.grafX(z*100, rho * 2e2, 'телеcкоп', 'z, mm', 'D, mm', legend)
else:
    legend = "F1", "F2"
    l = list(map(lambda l2,l3: (l2 * 100, (l2 + l3) * 100), l2[0], l3[0]))
    o.grafX(rho[:][-1] * 200, l,  'диаметр пучка от положения линз F2 и F3 от F1. F1 = F-74x2, F2 = 30, F3 = -74', 'D, mm', 'z, mm', legend)
