import matplotlib.pyplot as plt
from numpy import array, linspace, zeros

from Numeric.Esquemas_numéricos import (RK4, Crank_Nicolson, Euler,
                                        Euler_inverso, LeapFrog)
from Physics.Oscilators import OscilatorX
from Mathematics.EDOS import Cauchy_Problem
## Variables temporales ##

T = 20                           # Duración de la simulación en segundos
dt = 0.1                       # Paso de integráción en segundos
n = int(T/dt)                    # Numero de pasos 
t  = linspace(0,T,n)             # Vector de instantes separados dt



## Condiciones iniciales ##

U0 = array([0,1])

methods = [Euler,RK4, Crank_Nicolson, Euler_inverso, LeapFrog]
lista = ['Euler','RK4','Crank Nicolson','Euler inverso', 'LeapFrog']


for j in range (5):

    U = Cauchy_Problem(OscilatorX,t,U0,methods[j])
    print(U[len(t)-1,:])

    plt.title(f'Oscilador integrado con {lista[j]} con dt = {dt} s y T = {T} s')
    plt.xlabel("Tiempo (s)")
    plt.ylabel("X")
    plt.plot(t,U[:,0])

    #plt.savefig('Plots/Hito /'+lista[j]+ ' ' + str(dt)+'.png')

    plt.show()
