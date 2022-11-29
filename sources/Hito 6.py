import matplotlib.pyplot as plt
from numpy import array, linspace, zeros, size
from random import random
from Numeric.Esquemas_numéricos import (RK4, Crank_Nicolson, Euler,
                                        Euler_inverso, leapfrog)
from Numeric.ERK import ERK
from Physics.Restricted3BP import CR3BP, Lagrange_Points_Calculation
from Mathematics.EDOS import Cauchy_Problem

## Temporal variables ##

T = 80               # Integration duration [s]
dt = 0.001           # Integration step [s]
n = int(T/dt)                    # Number of steps  
t  = linspace(0,T,n)                # Time array

######### Circular Restricted Three Body Problem Resolution #########

mu = 3.0039e-7  # Earth - Sun
#mu = 1.2151e-2 # Earth - Moon
#mu = 7.1904e-4 # Jupiter - Sun
#mu = 2.8571e-4 # Saturn - Sun
#mu = 2.366e-4  # Saturn - Titan

def F(U,t):

    return CR3BP(U, t, mu)

######### Lagrange Points #########

NL = 5     # Number of Lagrange Points

U0 = zeros([NL,4])  # Assigning the initial values for the system resolution

U0[0,:] = array([0.8, 0.6, 0, 0])
U0[1,:] = array([0.8, -0.6, 0, 0])
U0[2,:] = array([-0.1, 0, 0, 0])
U0[3,:] = array([0.1, 0, 0, 0])
U0[4,:] = array([1.01, 0, 0, 0])

LagrangePoints = Lagrange_Points_Calculation(U0, NL, mu)

######### Orbits around Lagrange points #########

U0LP = zeros(4)

eps = 1e-2*random()
SeLG = 5

U0LP[0:2] = LagrangePoints[SeLG-1,:] + eps
U0LP[2:4] = eps

methods = [Euler,RK4, Crank_Nicolson, Euler_inverso, leapfrog, ERK]

for j in range (size(methods)):

    U_LP = Cauchy_Problem(F, t, U0LP, methods[j])

    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(U_LP[:,0], U_LP[:,1],'-',color = "r")
    ax1.plot(-mu, 0, 'o', color = "g")
    ax1.plot(1-mu, 0, 'o', color = "b")
    for i in range(NL):
        ax1.plot(LagrangePoints[i,0], LagrangePoints[i,1] , 'o', color = "k")

    ax2.plot(U_LP[:,0], U_LP[:,1],'-',color = "r")
    ax2.plot(LagrangePoints[SeLG - 1,0], LagrangePoints[SeLG - 1,1] , 'o', color = "k")

    ax1.set_xlim(-2,2)
    ax1.set_ylim(-2,2)
    ax1.set_title("Vista del sistema orbital")
    ax2.set_title("Vista del punto de Lagrange")
    ax2.set_xlim(LagrangePoints[SeLG - 1,0]-0.5,LagrangePoints[SeLG - 1,0]+0.5)
    ax2.set_ylim(LagrangePoints[SeLG - 1,1]-0.5,LagrangePoints[SeLG - 1,1]+0.5)
    fig.suptitle(f"Tierra-Sol - CR3BP ({methods[j].__name__}) - Órbita alrededor de L2 con t = {T}" )
    for ax in fig.get_axes():
        ax.set(xlabel='x', ylabel='y')
        ax.grid()

    plt.show()