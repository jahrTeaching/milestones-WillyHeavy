
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

print('Por favor, defina las condiciones iniciales del problema en el siguiente orden: X0, Y0, Vx0, Vy0')

X0 = float(input())
Y0 = float(input())
Vx0 = float(input())
Vy0 = float(input())

U0 = np.array([X0,Y0,Vx0,Vy0])

def Fs(X,Y,Vx,Vy):
    return np.array([Vx,Vy,-X/(X**2+Y**2)**(3/2),-Y/(X**2+Y**2)**(3/2)])

print('Sus condiciones iniciales son: ' + str(U0))

print('Defina el paso de tiempo para la integración')

dt = float(input())

print('Defina cuántos segundos quiere que dure la simulación')

t = float(input())
n = int(t/dt) 
Dt = np.linspace(0,t,n)

print('Por último, escoja el método de integración para la órbita escribiendo el número de la opción:')
print('1- Euler')
print('2- Crank')
print('3- Runge Kutta de orden 4')

choice = int(input())

if choice == 1:

##  MÉTODO DE EULER ##

# U(n+1) = U(n) + dt*F(n)



    print('Has elegido Euler')

    F0 = Fs(X0,Y0,Vx0,Vy0)
    F  = F0
    U  = U0

    for i in range(1,n+1):

         U = U + dt*F
         [X,Y,Vx,Vy] = U
         F = Fs(X,Y,Vx,Vy)

         Xplot = [X0] + [U[0]]
         Yplot = [Y0] + [U[1]]

         plt.title('Órbita con Euler con dt = '+ str(dt) + ' y ' + str(t) + ' segundos de integración')
         plt.xlabel("X")
         plt.ylabel("Y")
         plt.plot(Xplot,Yplot,'b.')


elif choice == 2:

    print('Has elegido Crank')

    F0 = Fs(X0,Y0,Vx0,Vy0)
    F  = F0
    U  = U0
    def sistema(Un1):

        x,y,z,s = Un1
        return(x - U[0] - dt/2*U[2] - dt/2*z,
               y - U[1] - dt/2*U[3] - dt/2*s, 
               z - U[2] -dt/2*(-U[0]/(U[0]**2 + U[1]**2)**(3/2)) - dt/2*(-x/(x**2 + y**2)**(3/2)),
               s - U[3] -dt/2*(-U[1]/(U[0]**2 + U[1]**2)**(3/2)) - dt/2*(-y/(x**2 + y**2)**(3/2)))
 

    for i in range(1,n+1):

        x,y,z,s = fsolve(sistema, (U[0],U[1],U[2],U[3]))

        U = np.array([x,y,z,s])

        Xplot = [X0] + [U[0]]
        Yplot = [Y0] + [U[1]]

        plt.title('Órbita con Crank-Nicolson con dt = '+ str(dt) + ' y ' + str(t) + ' segundos de integración')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.plot(Xplot,Yplot,'bo')
        




elif choice == 3:
    print('Has elegido Runge Kutta de orden 4')

    F0 = Fs(X0,Y0,Vx0,Vy0)
    F  = F0
    U  = U0

    for i in range(1,n+1):

        [X,Y,Vx,Vy] = U
        F  = Fs(X,Y,Vx,Vy)
        k1 = F
        k2 = Fs(X + k1[0]*dt/2,Y + k1[1]*dt/2,Vx + k1[2]*dt/2,Vy + k1[3]*dt/2)
        k3 = Fs(X + k2[0]*dt/2,Y + k2[1]*dt/2,Vx + k2[2]*dt/2,Vy + k2[3]*dt/2)
        k4 = Fs(X + k3[0]*dt,Y + k3[1]*dt,Vx + k3[2]*dt,Vy + k3[3]*dt)
        U  = U + (dt/6)*(k1+2*k2+2*k3+k4)

        Xplot = [X0] + [U[0]]
        Yplot = [Y0] + [U[1]]

        plt.title('Órbita con RK4 con dt = '+ str(dt) + ' y ' + str(t) + ' segundos de integración')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.plot(Xplot,Yplot,'bo')






print('La posición final y la velocidad final del cuerpo será =' + str(U))

plt.show()
