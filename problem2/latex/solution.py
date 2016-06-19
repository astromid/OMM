# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 14:43:10 2016

@author: Yan
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import sin

#функция граничных условий
def u_t0(x,y):
    return sin(2*np.pi*x) * sin(np.pi*y)

#TMA - tridiagonal matrix algorithm - метод прогонки
def TMA(A, B, C, F, k1, k2, mu1, mu2):
    N = A.size
    w = np.zeros(N+1)
    alpha = np.zeros(N+1)
    beta = np.zeros(N+1)
    
    #прямой ход прогонки
    alpha[1] = k1
    beta[1] = mu1

    for n in range(1, N):
        alpha[n+1] = B[n] / (C[n] - alpha[n]*A[n])
        beta[n+1] = (A[n]*beta[n] + F[n]) / (C[n] - alpha[n]*A[n])
    
    #обратный ход
    w[N] = (mu2 + beta[N]*k2) / (1 - alpha[N]*k2)
    
    for n in range(N-1, -1, -1):
        w[n] = alpha[n+1]*w[n+1] + beta[n+1]

    return w
    
    
#шаги по времени и координатам
N1 = 20
N2 = 40
S = 100

#рассматриваемый промежуток времени
T = 0.07

h1 = 1 / N1
h2 = 2 / N2
tau = T / S

#коэффициенты гамма_1 и гамма_2
g1 = tau / h1**2
g2 = tau / h2**2

#массив для искомого решения
#индексы: время(k) - x(i) - y(j)
u = np.zeros((S+1, N1+1, N2+1))
#массив для хранения значений на промежуточном слое k+1/2
u12 = np.zeros((S+1, N1+1, N2+1))

#заполняем граничные условия
for i in range(0, N1):
    for j in range(0, N2):
        u[0][i][j] = u_t0(i*h1, j*h2)

#основной цикл расчета
for k in range(S):
    #шаг k -> k + 1/2
    A = np.zeros(N1)
    B = np.zeros(N1)
    C = np.zeros(N1)
    F = np.zeros(N1)
    
    k1 = 0
    k2 = 0
    mu1 = 0
    mu2 = 0
    
    for j in range(1, N2):
        for i in range(1, N1):
            A[i] = g1 / 2
            B[i] = g1 / 2
            C[i] = 1 + g1
            F[i] = (g2 / 2) * (u[k][i][j-1] + u[k][i][j+1]) + (1 - g2) * u[k][i][j]

        w = TMA(A, B, C, F, k1, k2, mu1, mu2)
        
        for i in range(N1+1):
            u12[k][i][j] = w[i]

    #шаг k + 1/2 -> k + 1
    A = np.zeros(N2)
    B = np.zeros(N2)
    C = np.zeros(N2)
    F = np.zeros(N2)
    
    #значения k1, k2, mu1, mu2 не поменялись
    for i in range(1, N1):
        for j in range(1, N2):
            A[j] = g2 / 2
            B[j] = g2 / 2
            C[j] = 1 + g2
            F[j] = (g1 / 2) * (u12[k][i-1][j] + u12[k][i+1][j]) + (1 - g1) * u12[k][i][j]

        w = TMA(A, B, C, F, k1, k2, mu1, mu2)
        
        for j in range(N2+1):
            u[k+1][i][j] = w[j]

#визуализация
ax = Axes3D(plt.figure())
X = np.arange(0, 1 + h1, h1)
Y = np.arange(0, 2 + h2, h2)
X, Y = np.meshgrid(X, Y)
Z = np.transpose(np.array(u[S]))

ax.plot_surface(X, Y, Z, rstride=1, cstride=1, color='1')
ax.set_zlim3d(-1, 1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u(x, y, {time})'.format(time=str(T)))
ax.view_init(elev=30, azim=230)

plt.show()