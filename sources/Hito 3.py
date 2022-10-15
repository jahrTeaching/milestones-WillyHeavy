## HITO 3 ## 

from math import log10

import matplotlib.pyplot as plt
from numpy import array, linspace, zeros

from Numeric.Error import Richardson, Temporal_convergence_rate
from Numeric.Esquemas_numéricos import (RK4, Crank_Nicolson, Euler,
                                        Euler_inverso)
from Physics.Órbitas import Kepler

## Variables temporales ##

T = 20                           # Duración de la simulación en segundos
dt = 0.001                        # Paso de integráción en segundos
n = int(T/dt)                    # Numero de pasos 
t  = linspace(0,T,n)             # Vector de instantes separados dt

t1 = t
t2 = linspace(0,T,2*n)

## Condiciones iniciales ##

U0 = array([1,0,0,1])

log_N = zeros(n)

for i in range (n):

    log_N[i] = log10(n)

methods = [Euler,RK4, Crank_Nicolson, Euler_inverso]
lista = ['Euler','RK4','Crank Nicolson','Euler inverso']

for j in range(3): 

    Error = Richardson(Kepler, methods[j], t, U0)
    Error_norm = (Error[:,0]**2 + Error[:,1]**2)**(1/2)

    plt.figure(1)
    plt.subplot(211)
    plt.title(f'Error de {lista[j]} con dt = {dt} s y T = {T} s')
    plt.plot(t, Error[:,0],"r", label = "X position")
    plt.plot(t, Error[:,1],"b", label = "Y position")
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Error")
    plt.legend(loc = "lower left")
    plt.subplot(212)
    plt.plot(t, Error_norm, 'c', label = 'Module')
    plt.xlabel("t")
    plt.ylabel("Error")
    plt.legend(loc ='lower right')
    plt.grid()
    plt.savefig('Plots/Hito 3/ Error ' + lista[j]+ ' ' + str(dt)+'.png')
    #plt.show()
    plt.close('all')

    [log_E, log_N, log_E_lineal, log_N_lineal, order] = Temporal_convergence_rate(Kepler, methods[j], t, U0)
    
    
    plt.plot(log_N, log_E, "b", label = lista[j])
    plt.plot(log_N_lineal, log_E_lineal, "r", label = 'Linear regression')
    plt.legend(loc ='lower left')
    plt.xlabel("log(N)")
    plt.ylabel("log(U2-U1)")
    plt.title(f'{lista[j]} , order = {order}')
    plt.plot()
    plt.grid()
    plt.savefig('Plots/Hito 3/ Conv '+lista[j]+ ' ' + str(dt)+'.png')
    #plt.show() 
    plt.close('all')


'''
methods = [Euler,RK4, Crank_Nicolson, Euler_inverso]
lista = ['Euler','RK4','Crank Nicolson','Euler inverso']

for j in range (4):


    U = Cauchy_Problem(Kepler,t,U0,methods[j])
    print(U[len(t)-1,:])

    plt.title(f'Kepler integrado con {lista[j]} con dt = {dt} s y T = {T} s')
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.plot(U[:,0],U[:,1])

    plt.savefig('Plots/Hito 3/'+lista[j]+ ' ' + str(dt)+'.png')

    plt.show()

'''
