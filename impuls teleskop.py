#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 13:37:27 2018

@author: raptev

- Диаметр пучка на входе в манипулятор: 14.3 мм
- Полный угол расходимости на входе в манипулятор: 8.6 мрад
- На выходе манипулятора коллимированный пучок 8-15 мм
- Рабочая область: 0-200 мм от упора
- Шаг изменения диаметра: 0.5 мм 
- Цифровое обозначение диаметра: 8, 9, 10, 11, 12, 13, 14, 15 мм
- Длины волны: 532 и 1064 нм (одна общая шкала)
- Ошибка относительно выставленного диаметра: не более +/- 1 мм
- Диаметр пучка на оптических элементах должен быть не менее 7 мм
"""
from math import inf
import optics as o
import numpy as np
import time
def rangeZoom(D,f1,f2,f3):
    f11 = 1/(1/f1+1/f2)
    f12 = f3
    f21 = f1
    f22 = 1/(1/f2 + 1/f3)
    print(f11, f12, f21, f22)
    return D * f12 / f11, D * f22/ f21
start_time = time.clock()

n = 1, 2000# nx, ny
k = 2
f1 = -70e-2#74e-2
f2 = 58.2e-2#60e-2
f3 = -70e-2#4e-2
#f1 = -32.1e-2 #-74e-2/2
#f2 = 58.2e-2/2 #60e-2/2
#f3 = -32.1e-2#-74e-2/2
#f1 = -50E-3
#f2 = 30E-3
#f3 = 50e-3


#print('пучки', rangeZoom(14.3, f1,f2,f3))
y1 = o.OptY(*n, 14.3e-2/2, 8.6e-3)
l1 = 100e-2
l21 = 0 if (1 / (1 / f1 + 1 / f2) > 0) else f1 + f2
l22 = f1 + 1/(1/f2 + 1/f3)
L1 = o.OptSpace(*n,l1);

F1 = o.OptLens(*n, f1)
l2, y = o.vectorY(*n, l21, l22)
L2 = o.OptSpace(*n, l2)
F2 = o.OptLens(*n,f2)

opt = L1, F1, L2, F2
M = o.linMatrix(opt)
l3 = f3 + (M * y1).focus()
L3 = o.OptSpace(*n, l3)
F3 = o.OptLens(*n, f3)
opt2 = opt, L3, F3, L1

print (time.clock() - start_time, "seconds")
    
z, rho = o.lineOpt(opt2, y1)

D0 = 15e-2
i0=0
rho0 = 0

if (k == 0):
    legend = list(map(lambda y,l: 'L1 = ' + str(y*100) + 'mm' + ' L2 = ' + str(l * 100) + 'mm', y,l3[0]))
    o.grafX(z*100, rho * 2e2, 'телеcкоп', 'z, mm', 'D, mm', legend)
elif (k == 1):
    legend = "F2", "F3"
    l = list(map(lambda l2,l3: (l2 * 100, (l2 + l3) * 100), l2[0], l3[0]))
    o.grafX(rho[:][-1] * 200, l,  'диаметр пучка от положения линз F2 и F3 относительно F1.', 'D, mm', 'z, mm', legend)
elif (k == 2):
    for i, rhoi in enumerate(rho[-1][:]):
        if (abs(rhoi - D0/2) < .05e-2):
            i0 = i    
            rho0 = rhoi
            print('l1 = ', l2[0][i0] * 100, ' l2 = ', (l2[0][i0] + l3[0][i0]) * 100, ' D = ', rho0 * 200)
            break
print (time.clock() - start_time, "seconds")