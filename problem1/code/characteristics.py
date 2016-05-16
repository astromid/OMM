# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 19:43:00 2016

@author: Yan
"""
import matplotlib.pylab as plt
import numpy as np
from math import exp, pi, atan

def C(x):                                   #чтобы не считать каждый раз
    return (4/pi)*atan(x+2)

dx = 0.05                                   #шаг изменения параметра x0
dt = 0.05                                   #шаг изменения параметра t0
x0 = -1                                     #стартовое значение x0
t0 = 0                                      #стартовое значение t0
t = np.arange(0, 5, 0.05)                   #рассматриваемый отрезок времени и его дискретизация
plt.xlim(-1, 0)                             #настройки графика
plt.xlabel('x')
plt.ylabel('t')
for i in range(0, 20):                      #параметр t0 изменяется от 0 до 1
    x1 = []
    t0 = t0 + dt
    for curr_t in t:
        x1.append( (C(0)-2)*(curr_t - t0)*exp(-t0) )
    plt.plot(x1, t, 'r')
print('Final t0 = ', t0)   
for i in range(0, 20):                      #параметр x0 изменяется от -1 до 0
    x2 = []
    x0 = x0 + dx
    for curr_t in t:
        x2.append( x0 - curr_t*(2-C(x0)) )
    plt.plot(x2, t, 'b')
print('Final x0 = ',x0)
plt.show()