import matplotlib.pyplot as plt
from numpy import array, linspace, zeros, size

from Numeric.Esquemas_numéricos import (RK4, Crank_Nicolson, Euler,
                                        Euler_inverso, leapfrog)
from Physics.Oscilators import OscilatorX
from Mathematics.EDOS import Cauchy_Problem
from Numeric.Error import Characteristic_Polynomia

## Temporal variables ##

T = 20                           # Integration duration [s]
dt = 0.001                       # Integration step [s]
n = int(T/dt)                    # Number of steps  
t  = linspace(0,T,n)             # Time array

## Initial conditions for the oscilator ##

U0 = array([0,1])

methods = [Euler,RK4, Crank_Nicolson, Euler_inverso, leapfrog]
lista = ['Euler','RK4','Crank Nicolson','Euler inverso', 'LeapFrog']
figuras = ['euler','rk4','cranknicolson','Eulerinverso','leapfrog']
markers = ['r--','b--','g--']

for j in range (size(methods)):

    T = 20                                              # Integration duration [s]
    dt = array([0.1, 0.01, 0.001])                      # Integration step [s]
    

    for i in range(size(dt)):

        n = int(T/dt[i])                    # Number of steps  
        t  = linspace(0,T,n)                # Time array

        ### FIRST PART: OSCILATOR INTEGRATION ###

        U = Cauchy_Problem(OscilatorX,t,U0,methods[j])
        print(U[len(t)-1,:])

        plt.title(f'Oscilador integrado con {lista[j]}')
        plt.xlabel("Tiempo (s)")
        plt.ylabel("X",rotation = 0)
        plt.grid()
        plt.plot(t,U[:,0], markers[i], label = 'dt =' + str(dt[i]) + ' s')
        plt.legend(loc ='lower left')    
        
    plt.savefig('Plots/Hito 4/ osc' + figuras[j] +'.png')
    plt.close()
    
"""
    ### SECOND PART: SCHEME'S ABSOLUTE STABILITY REGION ###
    
    [A,B] = Characteristic_Polynomia(methods[j])
    X = linspace(-3,3,10)
    Y = linspace(-3,3,10)
    zer = zeros(10)
    fig = plt.figure()
    ax = fig.add_subplot()
    plt.title(f'Región de estabilidad absoluta de {lista[j]} ')
    plt.plot(X,zer,'k-')                         ## X axis
    plt.plot(zer,Y,'k-')                         ## Y axis
    plt.grid()
    plt.plot(A,B,'r',linewidth=4)                ## Stability region

    if methods[j] != Crank_Nicolson:
        ax.set_aspect('equal', adjustable = 'box')

    plt.xlabel("Re")
    plt.ylabel("Im",rotation = 0)
    plt.savefig('Plots/Hito 4/ ASR' + lista[j] + '.png')

"""