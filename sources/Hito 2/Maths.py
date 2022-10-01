from operator import matmul
from numpy import array, zeros
from numpy.linalg import inv, norm

## OPERACIONES MATEMÁTICAS ##

 #  
 #
 #  Inputs: 
 #         U : Vector estado en tn
 #         dt: Paso de tiempo 
 #         F : Sistema del que se quiere obtener la matriz Jacobiana o que se quiere solucionar con Newton-Raphson
 #
 #  Return: 
 #
 #         J : Matriz jacobiana 
 #         newton : solución del sistema no linear por Newton-Raphson 

def jacobiano(F,U):

    dim = len(U)
    Dx = 1e-3
    jacobian = array(zeros((dim,dim)))

    for i in range(dim):

        xj = array(zeros(dim))
        xj[i] = Dx
        jacobian[:,i] = (F(U + xj) - F(U - xj))/(2 * Dx)
    
    return jacobian

def newton(F, U0):

    dim = len(U0)
    Dx = array(zeros(dim))
    b  = array(zeros(dim))
    U1 = U0
    
    eps = 1
    iteration = 0
    itmax = 10000

    while eps > 1e-8 and iteration <= itmax:
   
        J = jacobiano(F,U1)
        b = F(U1)
        Dx = matmul(inv(J),b)
        U  = U1 - Dx
        eps = norm(U - U1)
        U1  = U
        iteration = iteration + 1

        if iteration == itmax:
            print('Máximo número de iteraciones alcanzado')
    
    return U