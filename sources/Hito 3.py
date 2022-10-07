## HITO 3 ## 

from numpy import array, linspace
from Numeric.Esquemas_numéricos import Euler,RK4, Crank_Nicolson, Euler_inverso
from Physics.Órbitas import Kepler
from Mathematics.EDOS import Cauchy_Problem
import matplotlib.pyplot as plt

## Variables temporales ##

T = 20                           # Duración de la simulación en segundos
dt = 0.001                        # Paso de integráción en segundos
n = int(T/dt)                    # Numero de pasos 
t  = linspace(0,T,n)             # Vector de instantes separados dt

## Condiciones iniciales ##

U0 = array([1,0,0,1])
U = U0

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

