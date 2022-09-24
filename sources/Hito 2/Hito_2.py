
## HITO 2 ## 

from numpy import array, zeros, linspace
import matplotlib.pyplot as plt


## Variables temporales ##

dt = 0.01
t  = 200
n  = int(t/dt) 
Dt = linspace(0,t,n)

## Condiciones iniciales ##

U0 = array([0,0,0,0])
U  = U0

def Euler(F,t,dt):

  U  = U0
  x  = array(zeros(n))
  y  = array(zeros(n))

  for i in range(1,n+1):

      U  = U + dt*F

      x[i] = U[0]
      y[i] = U[1]

      return U 

  plt.plot(x,y,'bo')


def Kepler(U,t):

    x = U[0]; y = U[1]; dxdt = U[2]; dydt = U[3]
    d = (x**2 + y**2)**1.5

    return array([ dxdt, dydt, -x/d, -y/d])

Euler(Kepler(U,t),t,dt)

plt.show