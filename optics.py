#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 13:37:27 2018

@author: raptev
"""

from math import sqrt, inf, pi, nan
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator


HeNe = 632e-9
class OptQ:
    optType = 1
    def __repr__(self): # вывод названия класса
        return "optQ()"
    def __str__(self): # печать с помощью print()
        return  ('-' + 'rho = ' + str(self.rho * 1000) +
              'мм,\n R = ' + str(self.R * 1000) + 
              'мм, \n lambda = ' + str(self.lambdaP * 1e9) + 'нм')
    #инициализация, в явном виде не нунжна, но возможно надо будет подправить     
    def __init__(self, rho = 1,R = inf, lambdaP = HeNe): 
        self.lambdaP = lambdaP
        self.rho = rho
        self.R = R
        try:
            self.q = 1. /(1./R + 1j * self.lambdaP/(pi * rho**2))
        except ZeroDivisionError:
            #print('! rho = ', rho, 'R = ', R)
            self.q = 0
        #self.focus()
#        self.focusD = 0
#        self.focusRho = 0
#        self.focusQ = 0
    def updateQ(self): # вычисление параметров пучка из комплекстных
        q = self
#        print(q)
        nx, ny = np.shape(q.q)
        q.rho = np.zeros((nx,ny))
        rho2 = np.zeros((nx,ny))
        q.R = np.zeros((nx,ny))
        for i in range(nx):
            for j in range(ny):
                try:
                    rho2[i][j] = q.lambdaP/ ((1./q.q[i][j]).imag * pi)
                    if rho2[i][j] >= 0 and rho2[i][j] < inf:
                        q.rho[i][j] = np.sqrt(rho2[i][j])
                    else:
                        q.rho[i][j]=0;
                   
                    q.R[i][j] = 1 / ((1/q.q[i][j]).real)
                except ZeroDivisionError:
                    q.q[i][j] = 0
                    q.rho[i][j] = 0
                    q.R[i][j] = 0
    #вычисление положение фокуса по параметрам пучка, не работает    
    def focus(self):
        try:
            self.focusD = -pi**2 * self.rho**4 * self.R / (pi**2 * self.rho**4 + 
                                             self.lambdaP**2 * self.R**2)
            q2 =  self.q + self.focusD
            rho2 = self.lambdaP / ((1./q2).imag * pi)
            if rho2 >= 0:
                rho = sqrt(rho2)
            else:
                rho=0;
        except ZeroDivisionError:
            self.focusD = 0
            rho = 0
            q2 = self.q
        self.focusRho = rho
        self.focusQ = q2
        
    def rho0(d, rho2, lambdaP = HeNe): # диаметр переятки по пятну пучка
        a = pi**2 * rho2**4 - 4 * d**2 * lambdaP**2
        rho0 = [0] * 2
        if a>=0:
            b = pi * rho2**2
            c = sqrt(2 * pi)
            d1 = sqrt(a) + b
            d2 = b - sqrt(a)
            if d1>=0:
                rho0[0] = sqrt(d1)/c
            if d2>=0:
                rho0[1] = sqrt(d2)/c
            return rho0
     
     
class OptM:
    optType = 0
    def __repr__(self):
        return "OptM()"
    def __str__(self):
        return  ('-' + str(self.M))
    def __init__(self, A=0,B=0,C=0,D=0):
        self.M = np.array([[A,B],[C,D]],float)
    def setOpt(self, A=0,B=0,C=0,D=0):
        self.M = np.array([[A,B],[C,D]],float)
    def abcdOpt(oM):
        nx, ny, tmp, tmp = np.shape(oM.M)
        A = np.zeros((nx,ny))
        B = np.zeros((nx,ny))
        D = np.zeros((nx,ny))
        C = np.zeros((nx,ny))
        for i in range(nx):
            for j in range(ny):
                A[i][j] = oM.M[i][j][0][0]
                B[i][j] = oM.M[i][j][0][1]
                C[i][j] = oM.M[i][j][1][0]
                D[i][j] = oM.M[i][j][1][1]
        return A, B, C, D
    def __mul__(self,M2):
        #перегрузка по типу
        if M2.optType == 0:        #произведение матриц
            nx, ny, tmp, tmp = np.shape(self.M)
            nx2, ny2, tmp, tmp = np.shape(M2.M)
            if nx != nx2 and ny != ny2:
               raise IOError("не верный размер массива!")
            M3 = OptM();
            M3.M = [0] * nx
            for i in range(nx):
                M3.M[i] = [0] * ny
                for j in range(ny):
                    M3.M[i][j] = np.dot(self.M[i][j], M2.M[i][j])
            return M3
        elif M2.optType == 1:   #вычисление пучка на выходе системы
            q1 = M2
            q2 = OptQ(lambdaP = q1.lambdaP)
            nx, ny, tmp, tmp = np.shape(self.M)
            nx2, ny2 = np.shape(q1.q)
            if nx != nx2 and ny != ny2:
                raise IOError("не верный размер массива!")
            A, B, C, D = OptM.abcdOpt(self)
            q2.q = (A * q1.q + B) / (C * q1.q + D)
            q2.updateQ()
            return q2
    def inverseQ(self):
        M = OptM();
        M.M = np.linalg.inv(self.M)
        return M
    # вычисление параметров пучка на первом зеркале резонатора
    def resonatorMode(self,lambdaP = 635e-9):
        nx, ny, tmp, tmp = np.shape(self.M)
        A, B, C, D = OptM.abcdOpt(self)
        det = (1 - ((D + A) * (1. / 2.)) ** 2)
        rho = np.zeros((nx,ny))
        rho2 = np.zeros((nx,ny))
        R = np.zeros((nx,ny))
        def fRhoR(A,B,C,D,det):
#            print(A,det)
            if det<=0:
                rho = 0
                R = 0
            else:
                rho2 = lambdaP * B / (pi * sqrt(det))
                rho = sqrt(rho2)
                R = -2 * B / (D - A)   
            return rho, R
        for i in range(nx):
            for j in range(ny):
                if (det[i][j] <= 0):
                    rho[i][j] = 0
                    R[i][j] = 0
                else:
                    rho2[i][j] = lambdaP * B[i][j] / (pi * sqrt(det[i][j]))
                    rho[i][j] = sqrt(rho2[i][j])
                    R[i][j] = -2 * B[i][j] / (D[i][j] - A[i][j])    
        q = OptQ(rho, R, lambdaP)
        return q
# на случай, если есписок оптических элементов состоит из списков
def linResonator(res, linRes=[]):
    if (type(res) is tuple):
        for i in res:
            linResonator(i,linRes)
        return linRes
    else:
        linRes.append(res)
        return linRes

#обход резонатора    
def resonatorMatrix(res):
    linRes = linResonator(res,[])
    M = linRes[0]
    for i in range(1,len(linRes)):
        M = M * linRes[i]
    for i in range(len(linRes) - 2, 0, -1):
        M = M * linRes[i]
    return M

# однократный обод системы
def linMatrix(res):
    linRes = linResonator(res,[])
    M = linRes[-1]
    for i in range(len(linRes) - 2, -1, -1):
        M = M * linRes[i]
        #linRes[i].printM()
    return M


def vectorXY(nx, ny, x):
    vec = np.ones((nx, ny)) * x
    return vec

        
def vectorX(nx, ny, x1, x2, dim = 1):
    vec = np.zeros((nx,ny)) 
    x = np.linspace(x1 * dim, x2 * dim, nx)
    for i in range(nx):
        for j in range(ny):
            vec[i][j] = x[i]
    return vec, x
        
def vectorY(nx, ny, y1, y2, dim = 1):
    vec = np.zeros((nx,ny)) 
    y = np.linspace(y1 * dim, y2 * dim, ny)
    for i in range(nx):
        for j in range(ny):
            vec[i][j] = y[j]
    return vec, y
     

class OptSpace(OptM):
    def __init__(self, nx, ny, d,):
        #проверка на целое
        if type(d) is float or type(d) is int:
            d = vectorXY(nx,ny,d)
        #размерность массива
        #заполняю массив M
        self.M = [0] * nx
        for i in range(nx):
            self.M[i] = [0] * ny
            for j in range(ny):
                self.M[i][j] = np.array([[1,d[i][j]],[0,1]],float)
        self.d = d
        
        
class OptLens(OptM):
    def __init__(self, nx, ny, f):
        if type(f) is float or type(f) is int:
            f = vectorXY(nx,ny,f)
        self.M = [0] * nx
        for i in range(nx):
            self.M[i] = [0] * ny
            for j in range(ny):
                self.M[i][j] = np.array([[1,0],[-1./f[i][j],1]],float)
        self.f = f
        
class OptMirror(OptM):
    def __init__(self, nx, ny,r=inf, n=1,):
        if type(r) is float or type(r) is int:
            r = vectorXY(nx, ny, r)
        if type(n) is float or type(n) is int:
            n = vectorXY(nx, ny, n)
        self.M = [0] * nx
        for i in range(nx):
            self.M[i] = [0] * ny
            for j in range(ny):
                self.M[i][j] = np.array([[1,0],[-2. * n[i][j] / r[i][j],1]],float)
        self.r = r
        self.n = n

class OptRefraction(OptM):
    def __init__(self, nx, ny, d, n=1):
        if type(d) is float or type(d) is int:
            d = vectorXY(nx,ny,d)
        if type(n) is float or type(n) is int:
            n = vectorXY(nx,ny,n)
        self.M = [0] * nx
        for i in range(nx):
            self.M[i] = [0] * ny
            for j in range(ny):
                self.M[i][j] = np.array([[1,d[i][j]/n[i][j]],[0,1]],float)
        self.d = d
        self.n = n
         
class OptSphere(OptM):
    def __init__(self, nx, ny, r, n1 = 1, n2 = 1):
        if type(r) is float or type(r) is int:
            r = vectorXY(nx,ny,r)
        if type(n1) is float or type(n1) is int:
            n1 = vectorXY(nx,ny,n1)
        if type(n2) is float or type(n2) is int:
            n1 = vectorXY(nx,ny,n2)
        self.M = [0] * nx
        for i in range(nx):
            self.M[i] = [0] * ny
            for j in range(ny):
                self.M[i][j] = np.array([[1,0],
                                   [(n2[i][j] - n1[i][j]) / r[i][j],1]],float)
        self.r = r
        self.n1 = n1
        self.n2 = n2

def grafX(x,z, title = 'graf',xlabel = 'x',ylabel = 'y', legend = 'graf1'):
    plt.plot (x, z)
    plt.grid(True)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(legend)
    
def graf2d(x, y, z, title, xlabel, ylabel, nbins = 15, zMin='min', zMax='max'):
    if zMin == 'min':
        zMin = z.min()
    if zMax == 'max':
        zMax = z.max()
                       
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    levels = MaxNLocator(nbins=nbins).tick_values(zMin, zMax)
#    levels = MaxNLocator(nbins=15).tick_values(z.min(), z.max())
    cmap = plt.get_cmap('PiYG')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    
    fig, (ax0) = plt.subplots(nrows=1)
    
    im = ax0.pcolormesh(x, y, z, cmap=cmap, norm=norm)
    fig.colorbar(im, ax=ax0)
    ax0.set_title(title)
    ax0.set_xlabel(xlabel)
    ax0.set_ylabel(ylabel)
    fig.tight_layout()
    plt.show()
    
    # функция изщет распределение пучка вдоль оси z
# дописать!!!!!
def lineOpt(res, q, nSpace = 100):
    # обязательно nx = 1
    nx,ny = np.shape(q.q)
    if nx != 1:
        raise IOError("nx должен равняться 1")
    linRes = linResonator(res,[])
    zi = np.zeros(ny)
    z = []
    rho = []
    R = []
    M = OptMirror(nx,ny)
    #a = "Hello" if foo() else "Goodbye"
    for i in linRes:
        if (type(i) is OptSpace):
            l=i.d
            dl = l/nSpace
            for j in range(nSpace):
                M = OptSpace(nx,ny,dl) * M
                q2 = M * q
                zi = zi + dl[0]              
                z.append([0]*ny)
                rho.append([0]*ny)
                R.append([0]*ny)
                for k in range(ny):
                    z[-1][k] = zi[k]
                    if q2.rho[0][k] > 0:
                        rho[-1][k] = q2.rho[0][k]
                        R[-1][k] = q2.R[0][k]
        else:
            M = i * M
            q2 = M * q
            z.append([0]*ny)
            rho.append([0]*ny)
            R.append([0]*ny)
            for k in range(ny):
                if q2.rho[0][k] > 0:
                    z[-1][k] = zi[k]
                    rho[-1][k] = q2.rho[0][k]
                    R[-1][k] = q2.R[0][k]
                    
    z = np.array(z)
    rho = np.array(rho)
    R = np.array(R)
    return z, rho, R
    
    #def grafRes(x,y,z):
#    P=range(1,60)
#    D = []
#    for i in P:
#        ri = resonator(i)
#        z = ri[0]
#        D.append(ri[1])
#    
#    z = np.array(z)
#    D = np.array(D)
#    P = np.array(P)
#    z, P = np.meshgrid(z, P)
#    fig = pylab.figure()
#    axes = Axes3D(fig)
#    axes.plot_surface(z, P, D, rstride=4, cstride=4, cmap = cm.jet)
#    axes.set_xlabel('z, m')
#    axes.set_ylabel('P, W')
#    axes.set_zlabel('D, m')
#    pylab.show()