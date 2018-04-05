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
import numpy as np

#constants and simulation parameters 
T=          0.1 #thermodynamical temperature 
gamma =     1.6 #transverse field strenght factor
kB =        1.0 #Boltzman constant
J =         1.0 #interaction constant
beta =      1.0/kB/T 
N =         16 #N quantum spins
M =         16 #additional dimension of spins
NT =        1000000 #liczba kroków czasowych, w których losowany jest jeden spin
snapNT =    NT/1000 #a value of magnetization is saved every 1000 steps of the simulation
Gsteps =    15 #number of evaluation points of gamma for the |<s>| = f(gamma) plot
Tsteps =    30 #number of times temperature is decreased during annealing

#initialization for first simulation
spins = fun.spinsInit(N,M)
J =     fun.JInit(N)

#keeping track of magnetization
globals.magnetizationInit()
globals.magnetizationHistoryInit()
globals.magnetization = spins.sum()

start = time.time();

'''Important!!!! Save time and remember to comment out the simulations. '''

#single simulation for arbitrary but constant parameters
fun.simulation(NT, snapNT, spins, J, N, M, beta, gamma)

plots.plot(spins, gamma, 10, T)
plots.plotMagHis(np.array(globals.magnetizationHistory)/(M*N))
print("magnetization", abs(globals.magnetization/(M*N)))
print("Avrage magnetization: ", fun.avrMag(globals.magnetizationHistory,N,M))



#simulation with decreasing gamma value, temperature is held constant
#gammaHistory = fun.simulationGammaRange(NT, snapNT, spins, J, N, M, beta, gamma, Gsteps)
#plots.plotGammaHis(gammaHistory)



#simulation of annealing (temperature varies, gamma is kept constant)
#fun.simulationAnnealing(NT, snapNT, spins, J, N, M, kB, gamma, T, Tsteps)

#simulation with varying N value

for m in [M * 2**i for i in range(5)]:
    nHistory = fun.simulationNRange(NT, snapNT, spins, J, m, beta, gamma, Gsteps)
    plots.plotNHistory(nHistory)


end = time.time();


print("time elapsed: ",end - start)





