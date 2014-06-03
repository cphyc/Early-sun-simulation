#!/usr/bin/python3
# -*- coding: utf-8 -*-
import numpy as np
from itertools import product
import matplotlib.pyplot as plt

verbose = False
debug = False
if debug:
    import ipdb

# constants 
Rsun = 696342e5           # cm
OmegaSun = 2.7e-6         # s⁻¹
Rgas = 8.314e7            # erg.K⁻¹.mol⁻¹
g = 27542.29              # cgs (taken at the surfce)

# set the parameters :
nu = 15                   # cm².s⁻¹
eta = 0 #496              # cm².s⁻¹
xi = 4e6                  # cm².s⁻¹
#xi = 4e3
r = 0.7*Rsun              # cm
theta = np.pi/4           # rad
R = r*np.cos(theta)       # cm
gamma = 5/3               # 
rho = 0.4                 # g.cm⁻³
T = 2.6e6                 # K
B_theta = 0               # Gauss (cgs)
v_a_theta = 0 #B_theta / np.sqrt(4*np.pi*rho)
                          # cm.s⁻¹

N2 = 6e-6
dPdr = -11239.4975692057  # dyn.cm⁻³
drho_gammaPdr = 1.121153213

def coeff (Omega, dlnOmegadlnr, k_R, k_Z):
    ''' Compute the coefficients αi as given in Menou et al. 2004 and 
    return them as an array.
    '''
    # recover dOmegadr w/ dlnOmegadlnr = r/Omega*dOmegadr
    dOmegadr = dlnOmegadlnr * Omega/r
    # precompute powers of k
    k2 = k_R**2 + k_Z**2
    k4 = k2**2
    k6 = k2**3
    k8 = k2**4
    k10 = k2**5

    # Some elements :

    # D(u) = kR/kZ du/dZ - du/dR
    #      = (kR/kZ cos(O) - sin(O)) * du/dr + (...) * du/dO
    # the right term is null for Omega (and we assume for P also)

    # u_r     = sin theta u_R + cos theta u_Z
    # u_theta = cos theta u_R - sin theta u_Z

    # precompute 1/gamma*rho * ...
    A = ( 1/(gamma*rho) * (k_R/k_Z/np.cos(theta) - 1/np.sin(theta))**2 * (-N2) )

    # precompute 
    #  1/R**3 D(R**4*Omega**2) = 2*R*Omega*(k_R/k_Z - 1)*dOmega/dr-4*Omega**2
    B = ( 2 *R*Omega* (k_R/k_Z/np.cos(theta) - 1/np.sin(theta)) * dOmegadr 
          - 4*Omega**2 )

    # precompute (k * v_a)**2 and its powers
    # k v_a = ktheta * Btheta / 4pirho
    #       = Btheta / ... * (kR cos theta - kZ sin theta)
    kva2 = v_a_theta * (k_R*np.cos(theta) - k_Z*np.sin(theta))
    kva4 = kva2**2

    # compute the coefficients
    a0 = k2/k_Z**2
    a1 = a0*(2*nu + 2*eta + xi)*k2
    a2 = ( a0*( k4*( nu**2 + eta**2 + 4*nu*eta + 
                     2*nu*xi + 2*eta*xi )
                + 2*kva2)
           - A - B )

    a3 = ( a0*( ( 2*eta*nu**2 + 2*nu*eta**2 + nu**2*xi + eta**2*xi 
                  + 4*nu*eta*xi ) * k6
                + 2*(nu + eta + xi)*k2 * kva2 )
           - (2*eta + nu)*k2 * A
           - (2*eta + xi)*k2 * B )

    a4 = ( a0*( (2*eta*xi*nu**2 + 2*nu*eta**2*xi + eta**2*nu**2)*k8
                + 2*(nu*eta + nu*xi + eta*xi)*k4*kva2
                + kva4)
           - ( (2*nu*eta*k4 + eta**2*k4 + kva2) * A )
           - ( (2*eta*xi*k4 + eta**2*k4 + kva2) * B )
           - 4*Omega**2*kva2 )

    a5 = ( a0*( xi*eta**2*nu**2*k10 + 2*xi*nu*eta*k6*kva2
                + xi*k2*kva4 )
           - ( (nu*eta**2*k6 + eta*k2*kva2) * A )
           - ( (xi*eta**2*k6 + xi*k2*kva2 ) * B )
           - 4*Omega**2*kva2*xi*k2 )

    if debug:
        ipdb.set_trace()

    return a0, a1, a2, a3, a4, a5

def loop(FGMs, FGMs_index, nOmega=10, ndOmega=10, nk=30, scale="log"):
    ''' Compute the Fastest Growing Modes (FGMs) over the Ω,∂lnΩ/∂lnr
    space by exploring the k_R, k_Z space.
    The FGM is at a given (Ω, ∂lnΩ/∂lnr):

    / σ  = -i*ω
    | αi = f(Ω, ∂lnΩ/∂lnr, k_R, k_Z)
    | α0 σ^5 + α1 σ^4 + α2 σ^3 + α3 σ^2 + α4 σ^1 + α5 = 0
    \ σFGM(Ω, ∂lnΩ/∂lnr) = max { σ_sol(k_R,k_Z) }
    '''
    Omega_range = [(31/(nOmega))*OmegaSun*x for x in range(nOmega)]
    dlnOmegadlnr_range = [-2.5/(ndOmega)*x for x in range(ndOmega)]
    # Omega_range = [(31*2/(nOmega))*OmegaSun*x for x in range(nOmega)]
    # dlnOmegadlnr_range = [-2.5*2/(ndOmega)*x for x in range(ndOmega)]

    # k ranges
    lmax = Rgas*T/g          # cm
    lmin = 1e5               # cm
    k_min, k_max = 2*np.pi/lmax, 2*np.pi/lmin
    
    if scale=="linear" or scale=="lin":
        # lin scale<<<<<<<<<<<<<<<<<<<<<
        k_range = [2*np.pi/lmax + n/(1.0*nk)*(2*np.pi/lmin-2*np.pi/lmax) 
                   for n in range(nk)]
        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    elif scale=="log":    
        # log scale ~~~~~~~~~~~~~~~~~~~~
        alpha = np.power(k_max/k_min, 1/(nk-1))
        k_range = [2*np.pi/lmax * alpha**n for n in range(nk)]
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    elif scale=="invert":
        # 1/l scale ^^^^^^^^^^^^^^^^^^^^
        l_range = [lmin + n/nk*(lmax-lmin) for n in range(nk)]
        k_range = [2*np.pi/l for l in l_range]
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    else : raise Exception("Invalid scale `%s`" % scale)

    # add the k < 0 values
    k_range += [-k for k in k_range]

    
    # iterate over the indexes of the Omega, dlnOmegadlnr ranges
    for a, b in product(range(nOmega), range(ndOmega)):
        FGM = 0
        Omega = Omega_range[a]
        dlnOmegadlnr = dlnOmegadlnr_range[b]
        max_k = (0,0)

        # iterate over all k-couples
        for c, d in product(range(nk), range(nk)):
            k_R = k_range[c]
            k_Z = k_range[d]

            p = coeff(Omega ,dlnOmegadlnr, k_R, k_Z)
            roots = np.roots(p)
            
            # local Fastest Growing Mode = biggest real value
            local_FGM = max(np.real(roots))

            # print(roots, local_FGM)
            
            if local_FGM > FGM:
                FGM = local_FGM
                max_k = (k_R, k_Z)
                if verbose:
                    print("# New FGM {}".format(FGM))

            if debug and Omega > 0 and dlnOmegadlnr < 0:
                ipdb.set_trace()
            if verbose:
                print(Omega,dlnOmegadlnr, k_R, k_Z, "\t", local_FGM)

        FGMs[a+nOmega*b] = [dlnOmegadlnr, Omega, FGM]
        FGMs_index[a,b] = FGM/OmegaSun

        if verbose:
            print("Ω: {:02f}, ∂lnΩ/∂lnr: {:02.2e}".format(Omega/OmegaSun, 
                                                          dlnOmegadlnr))

            print ("\t{:02.16e}\tk_R:".format(FGM)+
                   "{:02.2e}\tk_Z: {:02.2e}\n".format(max_k[0], max_k[1]))
    ext = (Omega_range[0]/OmegaSun, Omega_range[-1]/OmegaSun, 
           dlnOmegadlnr_range[0], dlnOmegadlnr_range[-1])
    return ext

