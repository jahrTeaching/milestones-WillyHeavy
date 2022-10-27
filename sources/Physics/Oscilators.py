from numpy import array
## OSCILATORS ##

 #  
 #
 #  Inputs: 
 #         U : Vector estado en tn
 #         t : tn 
 #
 #  Return: 
 #
 #         F(U,t) : Función derivada del vector de estado dU/dt = F(U,t)  

def OscilatorX(U,t):

    x    = U[0]
    dxdt = U[1]

    F = array([dxdt,-x])

    return F