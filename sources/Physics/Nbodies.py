from numpy import reshape, zeros, size
from numpy import linalg



def N_Bodies(U,t):  

    (Nb,Nc) = (2,3)    

    # Pinters
    Us = reshape(U,(Nb,Nc,2))           
    F  = zeros(len(U))
    Fs = reshape(F,(Nb,Nc,2))           # Respect Us order

    r  = reshape(Us[:,:,0],(Nb,Nc))     # Position, first index is the body and second its position
    v  = reshape(Us[:,:,1],(Nb,Nc))     # Velocity, first index is the body and second its velocity

    drdt = reshape(Fs[:,:,0], (Nb,Nc))
    dvdt = reshape(Fs[:,:,1], (Nb,Nc))

    dvdt[:,:] = 0                 

    for i in range(Nb):

        drdt[i,:] = v[i,:]

        for j in range(Nb):
            
            if j != i:                   # Only aplicable for different bodies

                d = r[j,:] - r[i,:]
                dvdt[i,:] += d[:]/linalg.norm(d)**3  

    return F