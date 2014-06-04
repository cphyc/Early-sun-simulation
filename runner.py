import numpy as np
from numpy import cos, sin
from itertools import product

verbose = False

# constants 
def get_param ():
    return lmin, lmax, OmegaSun

pi = 3.1415926535897932384626433832795028841971693
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
theta = pi/4              # rad
R = r*cos(theta)          # cm
gamma = 5/3               # 
rho = 0.4                 # g.cm⁻³
T = 2.6e6                 # K
B_theta = 0               # Gauss (cgs)
v_a_theta = 0 #B_theta / np.sqrt(4*pi*rho)
                          # cm.s⁻¹
N2 = 6e-6                 # s⁻²

# k ranges
lmax = Rgas*T/g           # cm
lmin = 1e5                # cm

# dPdr = -11239.4975692057  # dyn.cm⁻³
# drho_gammaPdr = 1.121153213

def  coeff (Omega, dlnOmegadlnr, k_R, k_Z):
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
    A = ( 1/(gamma*rho) * (k_R/k_Z/cos(theta) - 1/sin(theta))**2 * (-N2) )

    # precompute 
    #  1/R**3 D(R**4*Omega**2) = 2*R*Omega*(k_R/k_Z - 1)*dOmega/dr-4*Omega**2
    B = ( 2 *R*Omega* (k_R/k_Z/cos(theta) - 1/sin(theta)) * dOmegadr 
          - 4*Omega**2 )

    # precompute (k * v_a)**2 and its powers
    # k v_a = ktheta * Btheta / 4pirho
    #       = Btheta / ... * (kR cos theta - kZ sin theta)
    kva2 = v_a_theta * (k_R*cos(theta) - k_Z*sin(theta))
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

    ans = np.array([a0, a1, a2, a3, a4, a5])
    return ans

def  _loop(FGMs_index, Omega_range, dlnOmegadlnr_range, k_range,
			om_b, om_e,
			dom_b, dom_e):
    ''' Compute the Fastest Growing Modes (FGMs) over the Ω,∂lnΩ/∂lnr
    space by exploring the k_R, k_Z space.
    The FGM is at a given (Ω, ∂lnΩ/∂lnr):
    
    / σ  = -i*ω
    | αi = f(Ω, ∂lnΩ/∂lnr, k_R, k_Z)
    | α0 σ^5 + α1 σ^4 + α2 σ^3 + α3 σ^2 + α4 σ^1 + α5 = 0
    \ σFGM(Ω, ∂lnΩ/∂lnr) = max { σ_sol(k_R,k_Z) }
    '''
    nk = len(k_range)
    # iterate over the indexes of the Omega, dlnOmegadlnr ranges
    for a in range(om_b, om_e):
      for b in range(dom_b, dom_e):
        FGM = 0
        Omega = Omega_range[a]
        dlnOmegadlnr = dlnOmegadlnr_range[b]
        max_k_R = 0
        max_k_Z = 0

        # iterate over all k-couples
        for c in range(nk):
          for d in range(nk):
            k_R = k_range[c]
            k_Z = k_range[d]

            p = coeff(Omega ,dlnOmegadlnr, k_R, k_Z)
            roots = np.roots(p)
            
            # local Fastest Growing Mode = biggest real value
            local_FGM = max(np.real(roots))

            # print(roots, local_FGM)
            
            if local_FGM > FGM:
                FGM = local_FGM
                max_k_R = k_R
                max_k_Z = k_Z
                if verbose:
                    print("# New FGM {}".format(FGM))

            if verbose:
                print(Omega,dlnOmegadlnr, k_R, k_Z, "\t", local_FGM)

        FGMs_index[a,b] = FGM/OmegaSun

        if verbose:
            print("Ω: {:02f}, ∂lnΩ/∂lnr: {:02.2e}".format(Omega/OmegaSun, 
                                                          dlnOmegadlnr))

            print ("\t{:02.16e}\tk_R:".format(FGM)+
                   "{:02.2e}\tk_Z: {:02.2e}\n".format(max_k_R, max_k_Z))
    return FGMs_index

def loop(arg):
    # Unpack the parameters
    ( Omega_range, dlnOmegadlnr_range,
        k_range, om_b, om_e, dom_b, dom_e ) = arg
    
    # Create a receipter for the result
    FGMs_index = np.zeros((len(Omega_range),len(dlnOmegadlnr_range)))
    
    # Call the C _loop function
    return _loop(FGMs_index,
                 Omega_range,
                 dlnOmegadlnr_range,
                 k_range,                
                 om_b, om_e,
                 dom_b, dom_e)
