# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt

def plot(spins, gamma, sizeOfMatrix, T):
    fig=plt.figure()
    axis=fig.add_subplot(111)
    axis.imshow(spins, cmap='Greys', interpolation='none')#, origin='lower', extent=[gamma/sizeOfMatrix,gamma,T,T/sizeOfMatrix],aspect="auto")
    plt.title("moduł magnetyzacji")
    
    plt.show()

def plotMagHis(magnetizationHistory):
    plt.plot(magnetizationHistory)
    plt.ylabel('magnetization')
    plt.show()
    
def plotGammaHis(gammaHistory):
    plt.plot(gammaHistory[:,0], gammaHistory[:,1])
    plt.ylabel('|Average magnetization|')
    plt.xlabel('gamma')
    plt.show()
    
def plotGrid(array, gammaMin, gammaMax, TMin, TMax):
    fig, ax = plt.subplots()
    cax = ax.imshow(array, cmap='Greys', interpolation='none', origin='lower', extent=[gammaMin,gammaMax,TMin,TMax])
    plt.title("moduł magnetyzacji")
    plt.xlabel("T")
    plt.ylabel("Gamma")
    fig.colorbar(cax)
    plt.show()
    