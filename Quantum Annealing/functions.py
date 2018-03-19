# -*- coding: utf-8 -*-
import math
import numpy as np
import globals
import plots
#calculates the energy of the spin before the spin flip (nearest neighbour) 
def energyBeforeNN(spins, i, k, J, N, M, beta, gamma):
    energy = -spins[ i, k]*(J[ i, i - 1]*spins[i-1, k] + J[ i, (i + 1)%N] * spins[ (i+1)%N, k])/M - 0.5/beta * math.log(1.0/math.tanh(beta*gamma/M)) * spins[i,k]*(spins[i, k-1] + spins[i, (k+1)%M])
    return energy
 
#calculates the energy of the spin after the spin flip (nearest neighbour) 
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
    #energyAfter = energyArbitraryJ(spins, i, k, J, N, M, beta, gamma)
    
    deltaEnergy =  2 * energyAfter #- energyBefore
        
    
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

#calculates the average magnetiation of the system (averaging over the second half of the symylation) 
def avrMag(MagnetizationHistory, N, M):
    integral = 0
    for x in MagnetizationHistory[len(MagnetizationHistory)//2 : -1]:
        integral += x
    return integral/(N*M)/(len(MagnetizationHistory)//2)
    
#calculates the energy of the spin before the spin flip (works for a arbitrary matrix J_ij)
def energyArbitraryJ(spins, i, k, J, N, M, beta, gamma):
    temp = J[i,:]*spins[:,k]
    temp2 = temp.sum()
    energy = +spins[ i, k]*temp2/M + 0.5/beta * math.log(1.0/math.tanh(beta*gamma/M)) * spins[i,k]*(spins[i, k-1] + spins[i, (k+1)%M])
    return energy

def simulationGammaRange(NT, snapNT, spins, J, N, M, beta, gamma, Nsteps):
    gammaRange = np.linspace(0.01, gamma, Nsteps)
    gammaHistory = []
    for g in gammaRange:
        spins = spinsInit(N,M)
        J =     JInit(N)
        globals.magnetizationInit()
        globals.magnetizationHistoryInit()
        globals.magnetization = spins.sum()
        simulation(NT, snapNT, spins, J, N, M, beta, g)
        ar = np.array([g, abs(avrMag(globals.magnetizationHistory,N,M))])
        gammaHistory.append(ar)
    return np.array(gammaHistory)

def simulationAnnealing(NT, snapNT, spins, J, N, M, kB, gamma, T0, Tsteps):
    #initialization for first simulation
    spins = spinsInit(N,M)
    
    #keeping track of magnetization
    globals.magnetizationInit()
    globals.magnetizationHistoryInit()
    globals.magnetization = spins.sum()
    
    #calculating the decreas factor of temperature
    decreasFactor = (0.05/T0)**(1/Tsteps)
    
    #Array of slowly decreasing temperatures 
    iterationArray = np.zeros(Tsteps)
    for i in range(len(iterationArray)):
        iterationArray[i] = T0 * decreasFactor**i
    
    #Annealing
    for t in iterationArray:
        beta = 1/t/kB
        simulation(NT, snapNT, spins, J, N, M, beta, gamma)
        plots.plot(spins, gamma, 1, 1)

def simulationGammaRangeAnnealing(NT, snapNT, spins, J, N, M, kB, gammaMin, gammaMax, Nsteps, TMin, Tsteps, TMax, TIter):
    gammaRange = np.linspace(gammaMin, gammaMax, Nsteps)
    decreaseFactor = (0.05/TMin)**(1/Tsteps)
    iterationArray = np.zeros(Tsteps)   
    for i in range(len(iterationArray)):
        iterationArray[i] = TMin * decreaseFactor**(len(iterationArray)-i)    
    #print(iterationArray)
    grid = np.zeros((Nsteps, Nsteps), dtype=np.float)
    dT = (TMax-TMin)/TIter
    T = 0
    for g in range(len(gammaRange)):
        print("Gamma:", g, ", ", gammaRange[g])
        spins = spinsInit(N,M)
        J = JInit(N)
        globals.magnetizationInit()
        globals.magnetizationHistoryInit()
        globals.magnetization = spins.sum()
        #Annealing
        for T in range(TIter):
            #print("T0: ", TMin + T*dT)
            for t in iterationArray:
                beta = 1/t/kB
                simulation(NT, snapNT, spins, J, N, M, beta, gammaRange[g])
            grid[g,T] = abs(avrMag(globals.magnetizationHistory,N,M) )
            for t in iterationArray: 
                t += dT 
    #print(grid)
    return grid


