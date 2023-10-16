#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 15:29:19 2023

@author: alexandriaudenkwo
"""


from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir("/Users/alexandriaudenkwo/Project-1/phys-427-project-1")
result = np.loadtxt('output_blackhole_high_res.csv', delimiter=',')
print(len(result))
print(len(result[0]))
plt.imshow(result,cmap = 'inferno' ,origin="upper", extent=[-22,22,-11,11])
plt.xlabel('x/M', fontsize=12)
plt.ylabel('y/M', fontsize=12)

