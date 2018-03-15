# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 21:35:07 2018

@author: Paweł
"""

#libraries 
import functions as fun
import plots

#constants and simulation parameters 
T=      1.0 #termodynamic temperature 
gamma = 1.0 #transverse field strenght factor
k =     1.0 #Boltzman constant
J =     1.0 #interaction constant
beta =  1.0/k/T 
#matrix = 10 

N =     10 #N spinów kwantowych
M =     10 #additional dimension of spins
NT =    40000 #liczba kroków czasowych, w których losowany jest jeden spin


spins = fun.spinsInit(N,M)
J =     fun.JInit(N)

magnetization0 = spins.sum()


fun.simulation(NT, spins, J, N, M, beta, gamma)

magnetization = spins.sum()


plots.plot(spins, gamma, 10, T)



