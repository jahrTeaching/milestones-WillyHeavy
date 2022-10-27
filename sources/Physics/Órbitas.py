from numpy import array, reshape, zeros
from numpy import linalg

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


def N_Cuerpos(U,t):  

    (Nb,Nc) = (18,3)    # Tupla

    # Creación de punteros
    Us = reshape(U,(Nb,Nc,2))           # El 2 es posición y velocidad
    r  = reshape(Us[:,:,0],(Nb,Nc))     # Posición, el primer índice es el cuerpo y el segundo su posición
    v  = reshape(Us[:,:,1],(Nb,Nc))     # Velocidad, el primer índice es el cuerpo y el segundo su velocidad

    F  = zeros(len(U))
    Fs = reshape(F,(Nb,Nc,2))           # Tenemos que respetar el orden de la Us.
   
    drdt = reshape(Fs[:,:,0], (Nb,Nc))
    dvdt = reshape(Fs[:,:,1], (Nb,Nc))

    drdt[:,:] = v[:,:]                  #Esto rellena componentes de F porque apunta a F. 

    for i in range(Nb):

        dvdt[:,:] = 0

        for j in range(Nb):

            if j != i:

                d = r[j,:] - r[i,:]
                dvdt[i,:] += d[:]/linalg.norm(d)**3  #### Hay que hacer que i sea distinto de j

    return F


