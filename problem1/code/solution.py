# -*- coding: utf-8 -*-
"""
Created on Sun May 15 20:13:52 2016

@author: Yan
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import atan, exp, pi

#функции граничных условий
#на границе t = 0
def u_x0(x):
    return 2 - (4/pi)*atan(x+2)
    
#на границе x = 0
def u_0t(t):
    return (2 - (4/pi)*atan(2))*exp(-t)

#функция F из дивиргентной формы уравнения
def f_y(y):
    return -y**2 / 2
    
#производная f_y
def df_y(y):
    return -y

#генератор разностного уравнения на данном шаге
def gen_F(y_i1j, y_ij1, y_ij, h, tau, f_y=f_y):
    #разностное уравнение
    def F(y_i1j1):
        return (y_ij1 - y_ij + y_i1j1 - y_i1j) * h + \
            (f_y(y_i1j) - f_y(y_ij) + f_y(y_i1j1) - f_y(y_ij1)) * tau
    return F
    
#производная F
def dF(y, h, tau, df_y=df_y):
    return h + df_y(y) * tau
    
#итерационный метод поиска корня функции    
def newton_iteration(y0, F, dF, eps, h, tau):
    eps_cur = eps + 1
    while eps_cur > eps:
        y_iter = y0 - F(y0) / dF(y0, h, tau)
        eps_cur = np.abs(y_iter - y0)
        y0 = y_iter
        print("working in newton: eps_cur = " + str(eps_cur) + "esp = " + str(eps) + '\n')
    return y_iter
    
#шаги по времени и координате
N = 10
S = 10
T = 1

h = 1 / N
tau = T / S

#массив для искомого решения
y = np.zeros((S,N))

#заполняем граничные точки
for i in range(N):
    y[i][0] = u_x0(h*i)
    
for j in range(S):
    y[0][j] = u_0t(tau*j)

#основной цикл расчета
for j in range(S-1):
    for i in range(N-1):
        print(str(i) + ' ' + str(j) + " - i'm working...\n")
        F = gen_F(y[i+1][j], y[i][j+1], y[i][j], h, tau)
        y[i+1][j+1] = newton_iteration(y[i][j], F, dF, 0.0001, h, tau)
        
#визуализация
ax = Axes3D.Axes3D(plt.figure())
X = np.arange(-1, 0, h)
Y = np.arange(0, T, tau)
X, Y = np.meshgrid(X, Y)
Z = np.transpose(y)

ax.plot_surface(X, Y, Z, rstride=1, cstride=1, color='0.9')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u(x, t)')

plt.show()