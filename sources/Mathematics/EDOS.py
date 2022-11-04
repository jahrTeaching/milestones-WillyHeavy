from numpy import zeros, float64
from Numeric.Esquemas_numéricos import leapfrog


## Ecuaciones Diferenciales Ordinarias ##

 #  
 #
 #   Inputs: 
 #
 #         F(U,t) : Función derivada del vector de estado dU/dt = F(U,t)
 #         t : Vector de tiempo
 #         U : Vector estado en tn
 #         Esquema: Esquema temporal empleado para la integración 
 #
 #
 #
 #   Return: 
 #
 #         U Vector de estado en tn + dt  


def Cauchy_Problem(F,t,U0,Esquema):

    n, nv = len(t)-1,len(U0)

    U = zeros((n+1,nv), dtype=float64)

    dt = t[1] - t[0] 

    U[0,:] = U0

    if Esquema == leapfrog:

        U[1,:] = U[0,:] + dt*F(U[0,:],t[0])

        for i in range(1,n):

            U1 = U[i-1, :]
            U2 = U[i, :]
            U[i+1, :] = Esquema(U1, U2, dt, t[i], F)

    else: 
        for i in range(n):

            U[i+1,:] = Esquema(U[i,:], dt, t[i], F)

    return U


