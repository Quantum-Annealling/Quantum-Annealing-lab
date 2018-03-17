# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 21:35:07 2018

@author: Paweł
"""

#libraries 
import functions as fun
import plots
import time
import globals

#constants and simulation parameters 
T=      0.001 #termodynamic temperature 
gamma = 1.25 #transverse field strenght factor
k =     1.0 #Boltzman constant
J =     1.0 #interaction constant
beta =  1.0/k/T 
N =     15 #N spinów kwantowych
M =     15 #additional dimension of spins
NT =    100000 #liczba kroków czasowych, w których losowany jest jeden spin
snapNT =NT/1000
reductionFactot = 1.0/1.2 

spins = fun.spinsInit(N,M)
J =     fun.JInit(N)

#keeping track of magnetization
globals.magnetizationInit()
globals.magnetizationHistoryInit()
globals.magnetization = spins.sum()

start = time.time();
fun.simulation(NT, snapNT, spins, J, N, M, beta, gamma)
end = time.time();


plots.plot(spins, gamma, 10, T)
plots.plotMagHis(globals.magnetizationHistory)

print("time elapsed: ",end - start)

print("magnetization", abs(globals.magnetization/(M*N)))




