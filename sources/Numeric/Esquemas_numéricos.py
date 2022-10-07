import sys
sys.path.append('.')

from Mathematics.Maths import newton

## ESQUEMAS NUMÉRICOS ##

 #  
 #
 #  Inputs: 
 #         U : Vector estado en tn
 #         dt: Paso de tiempo
 #         t : tn 
 #         F(U,t) : Función derivada del vector de estado dU/dt = F(U,t)
 #
 #  Return: 
 #
 #         U Vector de estado en tn + dt  

        
def Euler(U,dt,t,F):

    return U + dt * F(U,t)

def RK4(U,dt,t,F):

    k1 = F(U,t)
    k2 = F(U + k1 * dt/2, t + dt/2)
    k3 = F(U + k2 * dt/2, t + dt/2)
    k4 = F(U + k3 * dt, t + dt)

    return U + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)

def Crank_Nicolson(U, dt, t, F):
    
    def Ecuación_Crank(Un1):
        return Un1 - U - dt/2 * (F(U,t) + F(Un1,t))

    return newton(Ecuación_Crank,U)

def Euler_inverso(U, dt, t, F):

    def Ec_Eulerinverso(Un1):
        return Un1 -U -dt*F(Un1,t)

    return newton(Ec_Eulerinverso,U)
