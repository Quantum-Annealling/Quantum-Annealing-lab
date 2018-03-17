# -*- coding: utf-8 -*-
import math
import numpy as np
import globals
#calculates the energy of the spin before the spin flip 
def energyBeforeNN(spins, i, k, J, N, M, beta, gamma):
    energy = -spins[ i, k]*(J[ i, i - 1]*spins[i-1, k] + J[ i, (i + 1)%N] * spins[ (i+1)%N, k])/M - 0.5/beta * math.log(1.0/math.tanh(beta*gamma/M)) * spins[i,k]*(spins[i, k-1] + spins[i, (k+1)%M])
    return energy
 
#calculates the energy of the spin after the spin flip 
def energyAfterNN(spins, i, k, J, N, M, beta, gamma):
    energy = +spins[ i, k]*(J[ i, i - 1]*spins[i-1, k] + J[ i, (i + 1)%N] * 
                       spins[ (i+1)%N, k])/M + 0.5/beta * math.log(1.0/math.tanh(beta*gamma/M)) * spins[i,k]*(spins[i, k-1] + spins[i, (k+1)%M])
    return energy

def flipSpin(spins, J, N, M, beta, gamma):
    
    #choose random spin on NxM grid
    i = np.random.randint(0, N, size = 1)
    k = np.random.randint(0, M, size = 1)
    
    #calculates the energy difference for the Metropolis algorithm
    
    #energyBefore = energyBeforeNN(spins, i, k, J, N, M, beta, gamma)
    energyAfter  = energyAfterNN(spins, i, k, J, N, M, beta, gamma)
    
    deltaEnergy =  2* energyAfter #- energyBefore
        
    
    #Metropolis algorithm
    if deltaEnergy < 0:
        spins[i,k] = -1 * spins[i,k]
        #keeping track of magnetization
        globals.magnetization = globals.magnetization + 2 * spins[i,k]
    else:
        boltzman = np.exp(-1.0*deltaEnergy*beta)
        uniform = np.random.random()
        if uniform < boltzman:
            #zmiana_magnetyzacji=-2*spins[x,y]
            #magnetyzacja+=zmiana_magnetyzacji
            spins[i,k]  = spins[i,k]*-1
            #keeping track of magnetization
            globals.magnetization = globals.magnetization + 2 * spins[i,k]

#initialization of grid of spins
def spinsInit(N, M):
    spins = np.random.randint(0,2,size=(N,M))*2-1 
    return spins

def JInit(N):
    J = np.zeros((N,N))
    for i in range(N):
        J[i, (i+1) % N] = 1
        J[(i+1) % N, i] = 1
    return J 
#function uses Monthe Carlo method for searching the minimum of the effective hamiltonian
def simulation(NT, snapNT, spins, J, N, M, beta, gamma):
    for x in range(NT):
        flipSpin(spins,J,N,M,beta,gamma)
        if x % snapNT  == 0:
            globals.magnetizationHistory.append(globals.magnetization)
    
    
    












