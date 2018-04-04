# -*- coding: utf-8 -*-
"""
Created on Wed Apr 04 14:48:51 2018
@author: Jakub
"""
import matplotlib.pyplot as plt
import numpy as np
import time

def readArray(input_file):
    t0 = time.time()
    f = open(input_file, 'r')
    lines = f.readlines()
    lines = [x.strip() for x in lines] 
    lines = [float(x) for x in lines]  
    f.close()  
    print('reading', str(input_file), "time= " , time.time() - t0)
    return np.asarray(lines)

magnetizationHistory = readArray('outputMagnetizationHistory.dat')
gammaHistory = readArray('outputGammaHistory.dat')
magnetizationGammaHistory = readArray('outputMagnetizationGammaHistory.dat')

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(7,6))
plt.locator_params(axis='y', nticks=6)
plt.locator_params(axis='x', nticks=10)
plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.7)
ax1.plot(magnetizationHistory, color='b', linestyle='-', label="magnetizationHistory")
ax1.legend(loc='upper right')
plt.grid(True)
plt.title('magnetization')
ax2.plot(gammaHistory, magnetizationGammaHistory, color='b', linestyle='-', label="gammaHistory")
ax2.legend(loc='upper right')
plt.grid(True)
plt.title('magnetization(gamma)')
plt.show()
