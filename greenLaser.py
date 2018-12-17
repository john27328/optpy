#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 13:37:27 2018

@author: raptev
"""
import gc
gc.collect()
from math import inf
import optics as o

rho1 = 3.57e-3
d1 = 1
#ищу радиус пучка в перетяжке по размеру пучка
rho0 = o.OptQ.rho0(d1, rho1, 1064e-9)
dl = 950e-3
f = 130e-3

q1 = o.OptQ(rho0[1],inf,532e-9)
q1.print('эксперимент')
M1 = o.OptSpace(dl)
M2 = o.OptLens(f)
q2 = M2 * M1 * q1       
q2.print('перетяжка')


q02 = o.OptQ(.16e-3/2, inf, 532e-9)
q02.print('теория1')
M1 = o.OptSpace(d1)
q2 = M1 * q02
q2.print('экран1')

M1 = o.OptSpace(dl)
M2 = o.OptLens(f)
q2 = M2 * M1 * q02      
q2.print('перетяжка1')
