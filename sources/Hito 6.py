import matplotlib.pyplot as plt

from numpy import array, linspace, zeros, size, around
from random import random
from Numeric.Esquemas_numéricos import (RK4, Crank_Nicolson, Euler,
                                        Euler_inverso, leapfrog)
from Numeric.ERK import ERK
from Physics.Restricted3BP import CR3BP, Lagrange_Points_Calculation, Lagrange_Points_Stability
from Mathematics.EDOS import Cauchy_Problem


## Temporal variables ##

T = 1000                                # Integration duration [s]
n = int(1e6)                            # Number of Points
t = linspace(0,T,n)                     # Time array


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
U0stabLP = zeros(4)
eps = 1e-3*random()
Lagrange_Points_List = array([1,2,3,4,5])
for k in range(5):

    selectedLP = k + 1

    if selectedLP == 5:
        label = 'L2'
    elif selectedLP == 4:
        label = 'L1'
    elif selectedLP == 3:
        label = 'L3'
    elif selectedLP == 2:
        label = 'L5'
    elif selectedLP == 1:
        label = 'L4'
    
    U0LP[0:2] = LagrangePoints[selectedLP-1,:] + eps
    U0LP[2:4] = eps

    U0stabLP[0:2] = LagrangePoints[selectedLP-1,:]
    U0stabLP[2:4] = 0

    eig = Lagrange_Points_Stability(U0stabLP, mu)
    print(around(eig.real,8))

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
        ax2.plot(LagrangePoints[selectedLP - 1,0], LagrangePoints[selectedLP - 1,1] , 'o', color = "k")

        ax1.set_xlim(-2,2)
        ax1.set_ylim(-2,2)
        ax1.set_title("Vista del sistema orbital")
        ax2.set_title("Vista del punto de Lagrange")
        ax2.set_xlim(LagrangePoints[selectedLP - 1,0]-0.02,LagrangePoints[selectedLP - 1,0]+0.02)
        ax2.set_ylim(LagrangePoints[selectedLP - 1,1]-0.02,LagrangePoints[selectedLP - 1,1]+0.02)
        fig.suptitle(f"Tierra-Sol - CR3BP ({methods[j].__name__}) - Órbita alrededor de {label} con t = {T}s" )
        for ax in fig.get_axes():
            ax.set(xlabel='x', ylabel='y')
            ax.grid()
            
        manager = plt.get_current_fig_manager()
        manager.full_screen_toggle()
        figure = plt.gcf()                  
        figure.set_size_inches(16, 8)       
        plt.savefig('Plots/Hito 6/ CR3BP ' + label +' '+ methods[j].__name__ +'.png', bbox_inches = 'tight')
        plt.close('all')
        #plt.show()

