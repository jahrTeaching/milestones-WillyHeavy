from numpy import array



## ÓRBITAS ##

 #  
 #
 #  Inputs: 
 #         U : Vector estado en tn
 #         t : tn 
 #
 #  Return: 
 #
 #         F(U,t) : Función derivada del vector de estado dU/dt = F(U,t)  




def Kepler(U,t):

    x = U[0]; y = U[1]; dxdt = U[2]; dydt = U[3]
    d = (x**2 + y**2)**1.5

    return array([ dxdt, dydt, -x/d, -y/d])