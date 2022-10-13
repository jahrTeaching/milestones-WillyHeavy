import sys

sys.path.append('.')

from numpy import linspace, log10, round_, zeros
from numpy.linalg import norm
from sklearn.linear_model import LinearRegression

from Mathematics.EDOS import Cauchy_Problem
from Numeric.Esquemas_numéricos import (RK4, Crank_Nicolson, Euler,
                                        Euler_inverso)


def Richardson(Problem,Scheme, t, U0):

    T  = t[-1]
    n  = len(t)
    t1 = t
    t2 = linspace(0,T,2*n)
    
    U1 = Cauchy_Problem(Problem, t1, U0, Scheme)
    U2 = Cauchy_Problem (Problem, t2, U0, Scheme )

    Error = zeros((n,len(U0)))
    
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


def Temporal_convergence_rate(Problem, Scheme, t, U0):  
    
    T  = t[-1]
    n  = len(t)
    t1 = t
    
    
    U1 = Cauchy_Problem(Problem, t1, U0, Scheme)
   

    m = 8
    log_E = zeros(m)
    log_N = zeros(m) 
    Error = zeros(m)
 
    for i in range (0,m):
 
        n = 2*n
        t2 = linspace(0,T,2*n)
        U2 = Cauchy_Problem (Problem, t2, U0, Scheme)

        Error[i] = norm(U2[int(n-1),:] - U1[int(n/2-1),:])
        log_E[i] = log10(Error[i])
        log_N[i] = log10(n)

        t1 = t2
        U1 = U2
        print(i)

    """   for j in range(0,m):

         if (abs(log_E[j]) > 12):

             break

    j = min(j, m)
    
    reg = LinearRegression().fit(log_N[0:j+1].reshape((-1, 1)),log_E[0:j+1])
    order = round_(abs(reg.coef_),1)

    log_N_lineal = log_N[0:j+1]
    log_E_lineal = reg.predict(log_N[0:j+1].reshape((-1, 1)))"""

    return [log_E, log_N]
   
   
    
    

