import matplotlib.pyplot as plt
from numpy import array, linspace, zeros, size

from Numeric.Esquemas_numéricos import (RK4, Crank_Nicolson, Euler,
                                        Euler_inverso, leapfrog)
from Physics.Oscilators import OscilatorX
from Mathematics.EDOS import Cauchy_Problem
from Numeric.Error import Characteristic_Polynomia
## Variables temporales ##

T = 20                           # Duración de la simulación en segundos
dt = 0.001                       # Paso de integráción en segundos
n = int(T/dt)                    # Numero de pasos 
t  = linspace(0,T,n)             # Vector de instantes separados dt



## Condiciones iniciales ##

U0 = array([0,1])


methods = [Euler,RK4, Crank_Nicolson, Euler_inverso, leapfrog]
lista = ['Euler','RK4','Crank Nicolson','Euler inverso', 'LeapFrog']

"""
for j in range (5):

    U = Cauchy_Problem(OscilatorX,t,U0,methods[j])
    print(U[len(t)-1,:])

    plt.title(f'Oscilador integrado con {lista[j]} con dt = {dt} s y T = {T} s')
    plt.xlabel("Tiempo (s)")
    plt.ylabel("X")
    plt.plot(t,U[:,0])

    #plt.savefig('Plots/Hito /'+lista[j]+ ' ' + str(dt)+'.png')

    plt.show()

"""

[A,B] = Characteristic_Polynomia(Euler_inverso)
maxA = max(abs(A))
maxB = max(abs(B))
X = linspace(-3,3,10)
Y = linspace(-3,3,10)
zer = zeros(10)
fig = plt.figure()
ax = fig.add_subplot()
plt.plot(X,zer,'k-')
plt.plot(zer,Y,'k-')
plt.grid()
plt.plot(A,B,'r',linewidth=4)
ax.set_aspect('equal', adjustable = 'box')
plt.show()
