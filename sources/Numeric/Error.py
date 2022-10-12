import sys
sys.path.append('.')

from numpy import linspace, zeros, log10, round_
from numpy.linalg import norm
from Mathematics.EDOS import Cauchy_Problem
from Numeric.Esquemas_numéricos import Euler, RK4, Crank_Nicolson, Euler_inverso
from sklearn.linear_model import LinearRegression


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
    t2 = linspace(0,T,2*n)
    
    U1 = Cauchy_Problem(Problem, t1, U0, Scheme)
    U2 = Cauchy_Problem (Problem, t2, U0, Scheme )
    
    m = 20 
    log_E = zeros(m)
    log_N = zeros(m) 
 
    for i in range (0,m):

        Error = norm(U2[2*i,:] - U1[i,:])
        log_E[i] = log10(Error)
        log_N[i] = log10(2*n)

    for j in range(0,n):

         if (abs(log_E[j]) > 12):

             break

    j = min(j, m)
    
    reg = LinearRegression().fit(log_N[0:j+1].reshape((-1, 1)),log_E[0:j+1])
    order = round_(abs(reg.coef_),1)

    log_N_lineal = log_N[0:j+1]
    log_E_lineal = reg.predict(log_N[0:j+1].reshape((-1, 1)))

    return [log_N, log_E, log_N_lineal, log_E_lineal]
   
   
    
    

