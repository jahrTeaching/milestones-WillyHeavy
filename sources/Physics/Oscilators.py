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

    return array(U[1],-U[0])