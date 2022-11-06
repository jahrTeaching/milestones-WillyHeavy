## HITO 5 ## 

import matplotlib.pyplot as plt
from numpy import array, linspace, zeros, size

from Numeric.Esquemas_numéricos import (RK4, Crank_Nicolson, Euler,
                                        Euler_inverso, leapfrog)
from Physics.Nbodies import N_Bodies
from Mathematics.EDOS import Cauchy_Problem

## Temporal variables ##

T = 20                           # Integration duration [s]
dt = array([0.1, 0.01, 0.001])   # Integration step [s]

## Initial conditions for the oscilator ##

U0 = array([1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0]) 

methods = [Euler,RK4, Crank_Nicolson, Euler_inverso, leapfrog]
lista = ['Euler','RK4','Crank Nicolson','Euler inverso', 'LeapFrog']
figuras = ['euler','rk4','cranknicolson','Eulerinverso','leapfrog']
markers = ['r--','b--','g--']

for j in range (size(methods)):

    T = 20                                              # Integration duration [s]
    dt = array([0.01, 0.001])                      # Integration step [s]

    for i in range(size(dt)):

        n = int(T/dt[i])                    # Number of steps  
        t  = linspace(0,T,n)                # Time array

        ### FIRST PART: OSCILATOR INTEGRATION ###

        U = Cauchy_Problem(N_Bodies,t,U0,methods[j])

        plt.plot(U[:,0], U[:,1], "b")
        plt.plot(U[:,6], U[:,7], "r")
        plt.show()

        fig = plt.figure()
        ax1 = fig.add_subplot(111,projection='3d')
        ax1.plot_wireframe(U[:,0].reshape((-1, 1)), U[:,1].reshape((-1, 1)), U[:,2].reshape((-1, 1)), color= "red")
        ax1.plot_wireframe(U[:,6].reshape((-1, 1)), U[:,7].reshape((-1, 1)), U[:,8].reshape((-1, 1)), color= "blue")
        plt.title(f'N cuerpos integrado con {lista[j]}')
        plt.xlabel("X")
        plt.ylabel("Y",rotation = 0)
        plt.grid()

        plt.show()
        plt.savefig('Plots/Hito 5/ N' + figuras[j] + str(dt[i])+'.png')
        

        
            
        
