import numpy as np
import matplotlib.pyplot as plt

print('Por favor, defina las condiciones iniciales del problema en el siguiente orden: X0, Y0, Vx0, Vy0')

X0 = float(input())
Y0 = float(input())
Vx0 = float(input())
Vy0 = float(input())

U0 = np.array([X0,Y0,Vx0,Vy0])

print('Sus condiciones iniciales son: ' + str(U0))
a = U0[0]+2
print(str(a))

print('Defina el paso de tiempo para la integración')

dt = float(input())

print('Defina cuántos segundos quiere que dure la simulación')

t = float(input())

n = int(t/dt) 

Dt = np.linspace(0,t,n)

##  MÉTODO DE EULER ##

# U(n+1) = U(n) + dt*F(n)

#U = np.array([X,Y,Vx,Vy])
#F = np.array([Vx,Vy,-X/(X^2+Y^2)^(3/2),-Y/(X^2+Y^2)^(3/2)])

F0 = np.array([Vx0,Vy0,-X0/(X0**2+Y0**2)**(3/2),-Y0/(X0**2+Y0**2)**(3/2)])
F  = F0
U  = U0

for i in range(1,n+1):

     U = U + dt*F
     [X,Y,Vx,Vy] = U
     F = np.array([Vx,Vy,-X/(X**2+Y**2)**(3/2),-Y/(X**2+Y**2)**(3/2)])
     Xplot = [X0] + [U[0]]
     Yplot = [Y0] + [U[1]]

     plt.plot(Xplot,Yplot,'bo')

print('La posición final y la velocidad final del cuerpo será =' + str(U))

plt.show()
