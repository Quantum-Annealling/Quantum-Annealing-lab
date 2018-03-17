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
    plt.show()