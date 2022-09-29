
## HITO 2 ## 

from numpy import array, zeros, linspace
from Esquemas_numéricos import Euler,RK4
import matplotlib.pyplot as plt


## Variables temporales ##

dt = 0.1
n  = 1000
t  = array (zeros(n))
T  = dt*n

## Condiciones iniciales ##

U = array([1,0,0,1])

x = array(zeros(n))
y = array(zeros(n))
x[0],y[0] = U[0],U[1]

def Kepler(U,t):

    x = U[0]; y = U[1]; dxdt = U[2]; dydt = U[3]
    d = (x**2 + y**2)**1.5

    return array([ dxdt, dydt, -x/d, -y/d])

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
plt.plot(x,y,'bo')

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
plt.plot(x,y,'bo')

plt.show()