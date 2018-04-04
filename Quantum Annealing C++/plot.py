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
plt.plot(magnetizationHistory)
plt.xlabel("iteracja")
plt.ylabel("magnetyzacja")
plt.show()
