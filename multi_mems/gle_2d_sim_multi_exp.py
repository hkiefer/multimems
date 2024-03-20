import numpy as np
import random
import math
import cmath 

from numba import njit


@njit()
def sqrtm(X):
    # Computing diagonalization
    evalues, evectors = np.linalg.eig(X)
    # Ensuring square root matrix exists
    assert (evalues >= 0).all()
    sqrt_matrix = evectors @ np.diag(np.sqrt(evalues)) @ np.linalg.inv(evectors)
    return np.real(sqrt_matrix)
    

@njit()
def integrate_gle_2_dim_exp(x0,v0, y0,ks, gammas, mass,nsteps, dt, kT,force_bins,force_matrix,fr_mode=1,uncorr_fr = False):

    n_dim = x0.shape[0]

    assert n_dim > 1, "n_dim must be larger than 1"

    n_exp = len(gammas)
    
    dot = np.dot
        
    x=np.zeros((int(nsteps),n_dim),dtype=np.float64)
  
    x[0]=x0
       
    inv_taus = np.zeros((gammas.shape))
    fac_rand = np.zeros((gammas.shape))
    
    for i in range(n_exp):
        
        inv_taus[i] = dot(np.linalg.inv(gammas[i]),ks[i])
        inv_gamma = np.linalg.inv(gammas[i])
        for k in range(n_dim):
            for l in range(n_dim):
                if inv_gamma[k][l] < 0:
                    inv_gamma[k][l] = 0
        if fr_mode == 1:
            
            fac_rand[i] = sqrtm(2*kT*np.linalg.inv(gammas[i])/dt)
        else:
            fac_rand[i] = np.sign(gammas[i])*np.sqrt(2*kT*inv_gamma/dt)
            

    if uncorr_fr:
        for l in range(n_exp):
            for i in range(0,n_dim):
                for j in range(0, n_dim):
                        if i != j:
                            fac_rand[l][i][j] = 0

                            
    xx=x0
    vv=v0
    yy=y0
    
    xi = np.zeros((n_exp,n_dim,))
    
    inv_mass = np.linalg.inv(mass)
    fr = np.zeros((n_exp,n_dim))
    
    for var in range(1,int(nsteps)):
        for i in range(n_exp):
            for j in range(n_dim):
                xi[i][j] = random.gauss(0.0,1.0) 

            fr[i] = dot(fac_rand[i],xi[i])

        kx1=dt*vv
        
        ff = dU_2D(xx,force_bins,force_matrix).reshape(n_dim,)
        
        kv1 = np.zeros((n_dim,))
        kv1-=dt*(dot(inv_mass,(ff)))
        
        ky1 = np.zeros((n_exp,n_dim,))
        y1 = np.zeros((n_exp,n_dim,))
        for i in range(n_exp):
            kv1 += dt*(-1*dot(inv_mass,dot(ks[i],(xx-yy[i]))))
            ky1[i] = dt*(-1*dot(inv_taus[i],(yy[i]-xx)) + fr[i])
            y1[i]=yy[i]+ky1[i]/2
        
        x1=xx+kx1/2
        v1=vv+kv1/2
        

        kx2=dt*v1
        
        ff = dU_2D(x1,force_bins,force_matrix).reshape(n_dim,)
        
        kv2 = np.zeros((n_dim,))
        kv2-=dt*(dot(inv_mass,(ff)))
        
        ky2 = np.zeros((n_exp,n_dim,))
        y2 = np.zeros((n_exp,n_dim,))
        for i in range(n_exp):
            kv2 += dt*(-1*dot(inv_mass,dot(ks[i],(x1-y1[i])))) 
            ky2[i] = dt*(-1*dot(inv_taus[i],(y1[i]-x1)) + fr[i])
            y2[i]=yy[i]+ky2[i]/2
        
        x2=xx+kx2/2
        v2=vv+kv2/2

        kx3=dt*v2
        
        ff = dU_2D(x2,force_bins,force_matrix).reshape(n_dim,)
        
        kv3 = np.zeros((n_dim,))
        kv3-=dt*(dot(inv_mass,(ff)))
        
        ky3 = np.zeros((n_exp,n_dim,))
        y3 = np.zeros((n_exp,n_dim,))
        for i in range(n_exp):
            kv3 += dt*(-1*dot(inv_mass,dot(ks[i],(x2-y2[i]))))
            ky3[i] = dt*(-1*dot(inv_taus[i],(y2[i]-x2)) + fr[i])
            y3[i]=yy[i]+ky3[i]
        
        x3=xx+kx3
        v3=vv+kv3

        kx4=dt*v3
       
    
        ff = dU_2D(x3,force_bins,force_matrix).reshape(n_dim,)
        
        kv4 = np.zeros((n_dim,))
        kv4-=dt*(dot(inv_mass,(ff)))
        
        ky4 = np.zeros((n_exp,n_dim,))
        y4 = np.zeros((n_exp,n_dim,))
        for i in range(n_exp):
            kv4 += dt*(-1*dot(inv_mass,dot(ks[i],(x3-y3[i]))))
            ky4[i] = dt*(-1*dot(inv_taus[i],(y3[i]-x3)) + fr[i])
            y4[i]=yy[i]+ky4[i]
            
            
        xx=xx + (kx1+2*kx2+2*kx3+kx4)/6
        vv=vv +(kv1+2*kv2+2*kv3+kv4)/6
        
        for i in range(n_exp):
            
            yy[i]+= (ky1[i]+(2*ky2[i])+(2*ky3[i])+ky4[i])/6
        

        x[var]=xx


    return x,vv,yy

@njit()
def integrate_langevin_multi_dim(x0,v0, gammas, mass,nsteps, dt, kT,force_bins,force_matrix,fr_mode=1,uncorr_fr = True):
    
    n_dim = x0.shape[0]

    assert n_dim > 1, "n_dim must be larger than 1"

    n_exp = len(gammas)
    
    dot = np.dot
        
    x=np.zeros((int(nsteps),n_dim),dtype=np.float64)
  
    x[0]=x0
       
    fac_rand = np.zeros((gammas.shape))
    
    for i in range(n_exp):
        
        gamma = gammas[i]
        for k in range(n_dim):
            for l in range(n_dim):
                if gamma[k][l] < 0:
                    gamma[k][l] = 0
        if fr_mode == 1:
          
            fac_rand[i] = sqrtm(2*kT*gammas[i]/dt)
        else:
            fac_rand[i] = np.sign(gammas[i])*np.sqrt(2*kT*gamma/dt)
            

    if uncorr_fr:
        for l in range(n_exp):
            for i in range(0,n_dim):
                for j in range(0, n_dim):
                        if i != j:
                            fac_rand[l][i][j] = 0

                            
    xx=x0
    vv=v0
    
    xi = np.zeros((n_exp,n_dim,))
    
    inv_mass = np.linalg.inv(mass)
    fr = np.zeros((n_exp,n_dim))
    
    for var in range(1,int(nsteps)):
        for i in range(n_exp):
            for j in range(n_dim):
                xi[i][j] = random.gauss(0.0,1.0) 

            fr[i] = dot(fac_rand[i],xi[i])

        kx1=dt*vv
        
        ff = dU_2D(xx,force_bins,force_matrix).reshape(n_dim,)
        
        kv1 = np.zeros((n_dim,))
        kv1-=dt*(dot(inv_mass,(ff)))
        
        for i in range(n_exp):
            kv1 += dt*(-1*dot(inv_mass,dot(gammas[i],(vv)) - fr[i]))
        
        x1=xx+kx1/2
        v1=vv+kv1/2
        

        kx2=dt*v1
        
        ff = dU_2D(x1,force_bins,force_matrix).reshape(n_dim,)
        
        kv2 = np.zeros((n_dim,))
        kv2-=dt*(dot(inv_mass,(ff)))
        
        for i in range(n_exp):
            kv2 += dt*(-1*dot(inv_mass,dot(gammas[i],(v1)) - fr[i]))
        
        x2=xx+kx2/2
        v2=vv+kv2/2

        kx3=dt*v2
        
        ff = dU_2D(x2,force_bins,force_matrix).reshape(n_dim,)
        
        kv3 = np.zeros((n_dim,))
        kv3-=dt*(dot(inv_mass,(ff)))
        
        for i in range(n_exp):
            kv3 += dt*(-1*dot(inv_mass,dot(gammas[i],(v2)) - fr[i]))
        
        x3=xx+kx3
        v3=vv+kv3

        kx4=dt*v3
       
    
        ff = dU_2D(x3,force_bins,force_matrix).reshape(n_dim,)
        
        kv4 = np.zeros((n_dim,))
        kv4-=dt*(dot(inv_mass,(ff)))
        
        for i in range(n_exp):
            kv4 += dt*(-1*dot(inv_mass,dot(gammas[i],(v3)) - fr[i]))
            
        xx=xx + (kx1+2*kx2+2*kx3+kx4)/6
        vv=vv +(kv1+2*kv2+2*kv3+kv4)/6
        
        
        x[var]=xx


    return x,vv


@njit()
def dU_2D(x,force_bins,force_matrix):

    value = np.zeros((len(x), 1), dtype=np.float64)

    #find indices once since they are the same for all force entries

    idxs = np.zeros((len(x)), dtype=np.int64)
    for j in range(len(x)):

        idx = bisection(force_bins[j],x[j])
        idxs[j] = int(idx)


    for i in range(len(x)):
   
        value[i] =  force_matrix[i][idxs[0]][idxs[1]]
       
    return value

@njit()
def bisection(array,value):
    '''Given an ``array`` , and given a ``value`` , returns an index j such that ``value`` is between array[j]
    and array[j+1]. ``array`` must be monotonic increasing. j=-1 or j=len(array) is returned
    to indicate that ``value`` is out of range below and above respectively.'''
    n = len(array)
    if (value < array[0]):
        return 0#-1
    elif (value > array[n-1]):
        return n-1
    jl = 0# Initialize lower
    ju = n-1# and upper limits.
    while (ju-jl > 1):# If we are not yet done,
        jm=(ju+jl) >> 1# compute a midpoint with a bitshift
        if (value >= array[jm]):
            jl=jm# and replace either the lower limit
        else:
            ju=jm# or the upper limit, as appropriate.
        # Repeat until the test condition is satisfied.
    if (value == array[0]):# edge cases at bottom
        return 0
    elif (value == array[n-1]):# and top
        return n-1
    else:
        return jl
    