# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 14:39:47 2018

@author: Jakub
"""



import numpy as np
import matplotlib.pyplot as plt
import math

k=1
J=1
T= 0.05
gamma = 0.27
beta=1/k/T
sizeOfMatrix = 15

N = 10 #NxN spinów w siatce
M = 10
NT=10000 #liczba kroków czasowych
snap_every=100

N_snapshots=NT//snap_every

spins= np.random.randint(0,2,size=(N,M))*2-1

magnetyzacja=spins.sum()

def energy(spins, beta, gamma, x, y, M):
    """
    Function receives a matrix describing states of the spins.
    The parameters beta and gamma are plain numbers with the 
    interpretation of temperature/nois in the system.
    """
    return  spins[x,y]*(spins[(x+1)%N, y]+spins[x-1, y])/M + 0.5/beta *math.log(1/math.tanh(beta*gamma/M)) * spins[x,y] * spins[x, (y+1)%M]
    

def plot(spins):
    fig=plt.figure()
    axis=fig.add_subplot(111)
    axis.imshow(spins, cmap='Greys', interpolation='none', origin='lower', extent=[gamma/sizeOfMatrix,gamma,T,T/sizeOfMatrix],aspect="auto")
    plt.title("moduł magnetyzacji")
    
    plt.show()
    

def flip_spin(spins, beta, gamma):
    """
    1. wybieramy losowo spin
        *losujemy x,y =random int[0,N]
    2.Liczymy enegie
        *E=-J *spin[x,y] * suma po sąsiadach z uwzglednieniem okresowych warunkow brzegowychs (PBC)
    3.liczymy praw. przerzucenia jako = exp(-2*E*beta)
    4. losujemy uniform(0,1)
    5. jesli boltzman > uniform(0,1): przerzucamy
    6.losowo wybrany spin jest przerzucony do gory nogami
    """
    global magnetyzacja
    x = np.random.randint(0, N, size = 1)
    y = np.random.randint(0, M, size = 1)
    #E=spins[x,y]*(spins[(x+1)%N, y]+spins[x-1, y]+spins[x, (y+1)%M]+spins[x, (y-1)%M])
    E = energy(spins, beta, gamma, x, y, M)
    
    boltzman=np.exp(-1*E*beta)
    uniform=np.random.random()
    if uniform<boltzman:
        zmiana_magnetyzacji=-2*spins[x,y]
        magnetyzacja+=zmiana_magnetyzacji
        spins[x,y]=spins[x,y]*-1
        

Beta = np.linspace(beta/sizeOfMatrix, beta, sizeOfMatrix)
Gamma =np.linspace(gamma/sizeOfMatrix, gamma, sizeOfMatrix)
Magnetyzacja = []
for k in Beta:
    m = []

    for j in Gamma:
        spins= np.random.randint(0,2,size=(N,M))*2-1
        magnetyzacja=spins.sum()
            
        #spins_history=np.empty((N_snapshots,N,M))
        mag_history=np.empty(NT)
        #plot(spins)
        for i in range(NT):
            #if(i%snap_every==0):
                #spins_history[i//snap_every]=spins
            mag_history[i]=magnetyzacja
            flip_spin(spins, k, j)
        m.append(float(magnetyzacja))
    Magnetyzacja.append(m)
            
def plot_mag():
    plt.plot(mag_history/(N*M))
    plt.show()
    

plot_mag()
plot(np.array(np.abs(Magnetyzacja)))

