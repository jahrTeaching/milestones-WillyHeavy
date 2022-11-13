import sys

sys.path.append('.')

from numpy import linspace, log10, round_, zeros, size, array, absolute, sqrt, float64
from numpy.linalg import norm
from sklearn.linear_model import LinearRegression
from cmath import  pi, sin, cos
from mpmath import findroot
from Mathematics.EDOS import Cauchy_Problem
from Numeric.Esquemas_numéricos import (RK4, Crank_Nicolson, Euler,
                                        Euler_inverso, leapfrog)

## Richardson Error calculation Function ##

def Richardson(Problem,Scheme, t, U0):

    n  = size(t)          # Número de puntos del mallado temporal
    T  = t[n-1]           # Tiempo final de la integración  
    t1 = t                # Primer mallado temporal
 
    t2 = linspace(0,T,2*n) # Segundo mallado temporal con el doble de puntos
    
    U1 = Cauchy_Problem(Problem, t1, U0, Scheme)
    U2 = Cauchy_Problem (Problem, t2, U0, Scheme )

    Error = zeros((n,size(U0)))
    
    if Scheme == Euler:
        q = 1
    elif Scheme == RK4:
        q = 4
    elif Scheme == Crank_Nicolson:
        q = 2
    elif Scheme == Euler_inverso:
        q = 1

    for i in range(n):
        Error[i,:] = (U2[2*i,:] - U1[i,:])/(1-1/(2**q))

    return Error 

## Schemes' convergence rate Function ##
def Temporal_convergence_rate(Problem, Scheme, t, U0):  
    
    
    n  = size(t)
    T  = t[n-1]
    t1 = t
    
    
    U1 = Cauchy_Problem(Problem, t1, U0, Scheme)
   
    m = 8               #Número de puntos que se quiere plotear.
    log_E = zeros(m)
    log_N = zeros(m) 
    Error = zeros(m)
    n = 2*n
    
 
    for i in range (0,m):

        t2 = linspace(0,T,(2**i)*n)
        U2 = Cauchy_Problem (Problem, t2, U0, Scheme)

        Error[i] = norm(U2[int((2**i)*n-1),:] - U1[int((2**i)*n/2-1),:])
        log_E[i] = log10(Error[i])
        log_N[i] = log10((2**i)*n)

        U1 = U2
        print(i)

        for j in range(0,m):

         if (abs(log_E[j]) > 12):

             break

    j = min(j, m)
    
    # Regresión lineal y cálculo de la pendiente 
    
    reg = LinearRegression().fit(log_N[0:j+1].reshape((-1, 1)),log_E[0:j+1]) 
    order = round_(abs(reg.coef_),1)

    log_N_lineal = log_N[0:j+1]
    log_E_lineal = reg.predict(log_N[0:j+1].reshape((-1, 1)))

    return [log_E, log_N, log_E_lineal, log_N_lineal, order]


## ABSOLUTE STABILITY REGION CALCULATION ##
def Stability_Region(Scheme): 

    N = 100
    x, y = linspace(-5, 5, 100), linspace(-5, 5, 100)
    rho =  zeros( (N, N),  dtype=float64)

    for i in range(N): 
      for j in range(N):

          w = complex(x[i], y[j])
          r = Scheme( 1., 1., 0., lambda u, t: w*u )
          rho[i, j] = abs(r) 

    return rho, x, y 


    
   
    
    

