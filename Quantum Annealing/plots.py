# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
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
    
def plotNHistory(nHistory):
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.set_xlabel('gamma/J')
    ax1.set_ylabel('|Magnetization|')

    #labels = []
    for i in range(np.shape(nHistory)[0]-1):
        ax1.plot(nHistory[0], nHistory[i+1], label= r'N = %i' % (16*2**i))
        #labels.append(r'$N = %i$' % (16*2**i))
    #plt.legend(labels, loc='upper center', bbox_to_anchor=[0.5, 1.1], columnspacing=1.0, labelspacing=0.0,handletextpad=0.0, handlelength=1.5,fancybox=True, shadow=False)
    ax1.legend(loc='upper right')
    plt.show()