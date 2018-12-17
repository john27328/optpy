#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 13:37:27 2018

@author: raptev
"""
from math import inf
import optics as o
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np

import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm

"""
    double d0=900e-6;//диаметр накачки, м
    f=1/.254/P*pow(d,2)/pow(d0,2); //qDebug()<<"f"<<f;
    linza(f);
)--l1--FNd--l2--l3--|
"""

rHr = -250e-3
l1 = 18e-3
P = 20
l2 = 25e-3
l3 = 185e-3

d0 = 900e-6
d = 1e-3

def resonator(P, l3):
    HR = o.OptMirror(rHr)
    L1 = o.OptSpace(l1)
    FNd = o.OptLens(1/.254/P * d**2/d0**2)
    L2 = o.OptSpace(l2)
    L3 = o.OptSpace(l3)
    OM = o.OptMirror()

    res = HR, L1, FNd, L2, L3, OM
    #res = HR, L1, FNd, L4, OM
    MNd = HR * L1 * FNd
    M1 = o.resonatorMatrix(res)
    q1 = M1.resonatorMode(532e-9)
    q2 = MNd * q1
#    q1.print()
#    q1.print()
#    z,rho,R = o.plotResonator(res,q1)
#    z = np.array(z)
#    R = np.array(R)
#    D = np.array(rho) * 2
    return q1, q2

    
def grafRes2d():
    nP = 500
    nL = 500
    P = np.linspace(1,100, nP)
    l3 = np.linspace(100, 600, nL) * 1e-3
    d = [0] * len(P)
    for i, Pi in enumerate(P):
        d[i] = [0] * len(l3)
        for j, li in enumerate(l3):
            q1, q2 = resonator(Pi,li)
            d[i][j]= q2.rho * 2
    
    x = np.array(l3)
    y = np.array(P)
    z = np.array(d)
    levels = MaxNLocator(nbins=15).tick_values(0, 1e-3)
#    levels = MaxNLocator(nbins=15).tick_values(z.min(), z.max())
    cmap = plt.get_cmap('PiYG')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    
    fig, (ax0) = plt.subplots(nrows=1)
    
    im = ax0.pcolormesh(x, y, z, cmap=cmap, norm=norm)
    fig.colorbar(im, ax=ax0)
    ax0.set_title('удлинение схемы magic3')
    ax0.set_xlabel('l3, m')
    ax0.set_ylabel('P, W')
    fig.tight_layout()
    plt.show()

    
grafRes2d()
