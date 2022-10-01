from numpy import zeros, float64


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

    U[0,:] = U0

    for i in range(n):

        U[i+1,:] = Esquema(U[i,:],t[i+1] - t[i],t[i],F)

    return U


