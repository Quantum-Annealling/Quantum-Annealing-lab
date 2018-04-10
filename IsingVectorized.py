from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt


def delta_energy(snew,sold,Jij):
    return np.dot(Jij,snew) - np.dot(Jij,sold)

def update_ising(si,Jij,beta=1.0):
    
    # vector of random {-1,1}
    snew = 2*np.random.randint(0,2,size=si.size)-1
    
    dE = delta_energy(snew,si,Jij)
    
    # metropolis update
    random_metropolis = np.random.uniform(0,1,size=si.size) # vector of random numbersfrom interval [0,1]
    ifupdate = np.zeros(si.size) 
    ifupdate[np.where(dE < 0)] = 1                 # for dE < 0 set all flips to one
    idx = np.where(dE > 0)[0]
    ifupdate[np.where( random_metropolis[idx] < np.exp(-beta*dE[idx]) )] = 1                 # for dE > 0 set all flips to one
    
    
    
    return si*(1-ifupdate)  +  snew*ifupdate       # for ifupdate == 0 -> no change
                                                   # for ifupdate == 1 -> change spin
def energy_ising(si,Jij):
    return np.dot( si, np.dot(Jij,si) )



if __name__ == '__main__':
    nt = 10000
    N = 1024
    J = 1.0
    beta = 3.0
    
    # construct interaction matrix
    Jij = np.zeros([N,N],dtype=np.float64)
    for i in range(N):
        Jij[i,i-1] = J; Jij[i,(i+1)%N] = J;
    
    s = 2*np.random.randint(0,2,size=N)-1
    
    for it in range(nt):
        s = update_ising(s,Jij,beta=beta)
        print(it, energy_ising(s,Jij), s.mean())
    
    plt.plot(s)
    plt.show()