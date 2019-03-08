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
l, y = o.vectorY(*n, 20* 1e-3, 80* 1e-3) 
L = o.OptSpace(*n,l);
F = o.OptLens(*n,30e-3)
y1 = o.OptY(*n, 8e-3, 0)
opt = L, F, L
#print(q1)
#print(q2)

#legend = list(map(lambda y: 'rHr = ' + str(y * 1000) + 'mm', y))
#
#o.grafX(x, q1.rho * 2e3, 'график 1', 'P, W', 'D, mm', legend)

print (time.clock() - start_time, "seconds")
    
#o.graf2d(x,y,q2.rho * 2e6, 'ris1', 'P, Вт', 'радиус зеркала, мм', 150, 200, 1000)
#print (time.clock() - start_time, "seconds")

z, rho= o.lineOpt(opt, y1)
legend = list(map(lambda y: 'P = ' + str(y) + 'W', y))
o.grafX(z, rho * 2e3, 'график 1', 'z, mm', 'D, mm',legend)
