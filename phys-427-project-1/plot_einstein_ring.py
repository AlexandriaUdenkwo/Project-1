#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 16:44:34 2023

@author: alexandriaudenkwo
"""

import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir("/Users/alexandriaudenkwo/phys-427-project-1")

result = np.loadtxt('output_einstein_ring_RKDP_4_pi.csv', delimiter=',')

plt.plot(result[:,1]*np.sin(result[:,2])*np.cos(result[:,3]), result[:,1]*np.sin(result[:,2])*np.sin(result[:,3]), label="Infinite well")

#Axis stuff
#gca means get current axes
ax = plt.gca()
ax.set_xlabel('x/M')
ax.set_ylabel('y/M')
plt.title("Einstein Ring")
plt.show()