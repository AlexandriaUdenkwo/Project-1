#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 22:50:31 2023

@author: alexandriaudenkwo
"""

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir("/Users/alexandriaudenkwo/phys-427-project-1")

resultA = np.loadtxt('output_case_A.csv', delimiter=',')
# resultB = np.loadtxt('output_case_B.csv', delimiter=',')
# resultC = np.loadtxt('output_case_C.csv', delimiter=',')
# resultD = np.loadtxt('output_case_D.csv', delimiter=',')
# resultE = np.loadtxt('output_case_E.csv', delimiter=',')


#plt.figure(1)
plt.plot(resultA[:,0],(np.abs(resultA[:,1] - (1+np.sqrt(2))))/(1+np.sqrt(2)) , label="error")


fig = plt.figure()
ax = plt.axes(projection='3d')




# zline = np.linspace(0, 15, 1000)
# xline = np.sin(zline)
# yline = np.cos(zline)
# ax.plot3D(xline, yline, zline, 'gray')

# plt.figure(2)
zline = resultA[:,1]*np.cos(resultA[:,2])
xline = np.sqrt(resultA[:,1]**2 + 1**2)*np.sin(resultA[:,2])*np.cos(resultA[:,3])
yline = np.sqrt(resultA[:,1]**2 + 1**2)*np.sin(resultA[:,2])*np.sin(resultA[:,3])
ax.plot3D(xline, yline, zline, 'blue')



plt.show()


# plt.figure(2)
# zline = resultB[:,1]*np.cos(resultB[:,2])
# xline = np.sqrt(resultB[:,1]**2 + 1**2)*np.sin(resultB[:,2])*np.cos(resultB[:,3])
# yline = np.sqrt(resultB[:,1]**2 + 1**2)*np.sin(resultB[:,2])*np.sin(resultB[:,3])
# ax.plot3D(xline, yline, zline, 'gray')



# zline = resultC[:,1]*np.cos(resultC[:,2])
# xline = np.sqrt(resultC[:,1]**2 + 1**2)*np.sin(resultC[:,2])*np.cos(resultC[:,3])
# yline = np.sqrt(resultC[:,1]**2 + 1**2)*np.sin(resultC[:,2])*np.sin(resultC[:,3])
# ax.plot3D(xline, yline, zline, 'gray')


# plt.figure(4)

# zline = resultD[:,1]*np.cos(resultD[:,2])
# xline = np.sqrt(resultD[:,1]**2 + 1**2)*np.sin(resultD[:,2])*np.cos(resultD[:,3])
# yline = np.sqrt(resultD[:,1]**2 + 1**2)*np.sin(resultD[:,2])*np.sin(resultD[:,3])
# ax.plot3D(xline, yline, zline, 'gray')

# plt.figure(5)
# zline = resultE[:,1]*np.cos(resultE[:,2])
# xline = np.sqrt(resultE[:,1]**2 + 1**2)*np.sin(resultE[:,2])*np.cos(resultE[:,3])
# yline = np.sqrt(resultE[:,1]**2 + 1**2)*np.sin(resultE[:,2])*np.sin(resultE[:,3])
# ax.plot3D(xline, yline, zline, 'gray')

plt.show()






