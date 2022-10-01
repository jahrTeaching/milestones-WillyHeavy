
## HITO 2 ## 

from numpy import array, zeros, linspace
from Esquemas_numéricos import Euler,RK4, Crank_Nicolson, Euler_inverso
from Órbitas import Kepler
from EDOS import Cauchy_Problem
import matplotlib.pyplot as plt

## Variables temporales ##

T = 100                          # Duración de la simulación en segundos
dt = 0.001                       # Paso de integráción en segundos
n = int(T/dt)                    # Numero de pasos 
t  = linspace(0,T,n)             # Vector de instantes separados dt

## Condiciones iniciales ##

U0 = array([1,0,0,1])
U = U0

x = array(zeros(n))
y = array(zeros(n))
x[0],y[0] = U[0],U[1]

## Euler ## 

for i in range(1,n):

    t[i] = dt*i
    U  = Euler(U,dt,t,Kepler)
    x[i] = U[0]
    y[i] = U[1]
       
print(U)

plt.title('Órbita con Euler con dt = '+ str(dt) + 's y ' + str(T) + ' segundos de integración')
plt.xlabel("X")
plt.ylabel("Y")
plt.plot(x,y)

plt.show()

U = array([1,0,0,1])

## Runge-Kutta 4 ## 

for i in range(1,n):

    t[i] = dt*i
    U  = RK4(U,dt,t,Kepler)
    x[i] = U[0]
    y[i] = U[1]

print(U)

plt.title('Órbita con RK4 con dt = '+ str(dt) + 's y ' + str(T) + ' segundos de integración')
plt.xlabel("X")
plt.ylabel("Y")
plt.plot(x,y)

plt.show()

U = array([1,0,0,1])

## Crank-Nicolson ## 

for i in range(1,n):

    t[i] = dt*i
    U  = Crank_Nicolson(U,dt,t,Kepler)
    x[i] = U[0]
    y[i] = U[1]
       
print(U)

plt.title('Órbita con CN con dt = '+ str(dt) + 's y ' + str(T) + ' segundos de integración')
plt.xlabel("X")
plt.ylabel("Y")
plt.plot(x,y)

plt.show()

U = array([1,0,0,1])

## Euler inverso ## 

for i in range(1,n):

    t[i] = dt*i
    U  = Euler_inverso(U,dt,t,Kepler)
    x[i] = U[0]
    y[i] = U[1]
       

plt.title('Órbita con Euler inverso con dt = '+ str(dt) + 's y ' + str(T) + ' segundos de integración')
plt.xlabel("X")
plt.ylabel("Y")
plt.plot(x,y)

print(U)

plt.show()




#def Cauchy_Problem(F,t,U0,Esquema):

#    n, nv = len(t)-1,len(U0)

#    U = zeros((n+1,nv), dtype=float64)

#    U[0,:] = U0

#    for i in range(n):

#        U[i+1,:] = Esquema(U[i,:],t[i+1] - t[i],t[i],F)

#    return U

#for j in range (2):


#    methods = [Euler,RK4]
#    U = Cauchy_Problem(Kepler,t,U0,methods[j])
#    print(U[len(t)-1,:])

#    plt.title('Solución al problema de Cauchy empleando '+ str(methods[j]) + 'con dt = '+ str(dt) + 's y ' + str(T) + ' segundos de integración')
#    plt.xlabel("X")
#    plt.ylabel("Y")
#    plt.plot(U[:,0],U[:,1])

#    plt.show()

##################
T = 100                          # Duración de la simulación en segundos
dt = 0.001                       # Paso de integráción en segundos
n = int(T/dt)                    # Numero de pasos 
t  = linspace(0,T,n)
U0 = array([1,0,0,1])

U = Cauchy_Problem(Kepler,t,U0,Euler_inverso)
print(U[len(t)-1,:])

plt.title('Solución al problema de Cauchy empleando  Euler con dt = '+ str(dt) + 's y ' + str(T) + ' segundos de integración')
plt.xlabel("X")
plt.ylabel("Y")
plt.plot(U[:,0],U[:,1])

plt.show()