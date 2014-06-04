import matplotlib.pyplot as plt
try:
    from crunner import loop, get_param
    print("Using Cython version")
except:
    from runner import loop, get_param
    print("Using Python version")
    
import numpy as np
try:
    from multiprocessing import Pool
    multiproc = 16
    pool = Pool(processes=multiproc)
    print ("Using %d processors" % multiproc)
except:
    multiproc = False
    
from itertools import product

# Extension of the omega domain and the k domain
nk = 10
nOmega = 10
ndOmega = nOmega

# logarithmic scale for k
lmin, lmax, OmegaSun = get_param()
l_range = [lmin + n/nk*(lmax-lmin) for n in range(nk)]
k_range = [2*np.pi/l for l in l_range]
k_range += [-k for k in k_range]

k_range = np.array(k_range)

# Omega ranges
Omega_range = np.array([31/nOmega*OmegaSun*i for i in range(nOmega)])
dlnOmegadlnr_range = np.array([-2.5/ndOmega*i for i in range(ndOmega)])

if multiproc:
    params = []
    for nom in range(nOmega//4 + 1):
        for ndom in range(ndOmega//4 + 1):
            om_b = 4*nom
            om_e = min(om_b+4, nOmega)
            dom_b= 4*ndom
            dom_e= min(dom_b+4,ndOmega)

            params.append( [Omega_range, dlnOmegadlnr_range,
                            k_range, om_b, om_e, dom_b, dom_e] )
    results = pool.map(loop, params)
else:
    results = [loop([Omega_range, dlnOmegadlnr_range,
                    k_range, 0, nOmega, 0, ndOmega])]

# # grep the input by blocks of 4 Omega and 4 dlnOmegadlnr
# gen1 = ( (Omega_range[4*n:min(4*n+4,nOmega)], 4*n )
#         for n in range(nOmega//4+1) )
# gen2 = ( (dlnOmegadln_range[4*n:min(4*n+4,ndOmega)], 4*n)
#         for n in range(ndOmega//4+1) )
# k_range = np.array(k_range)

# # creates a generator that grep all the 4*4 Omegas possible
# gen = ([FGMs, FGMs_index, om, dom, k_range, "log", om_offset, dom_offset]
#        for (om, om_offset), (dom, dom_offset) in product(gen1,gen2))

# # run the loop we get in return 
# results = pool.map(loop, gen)

# get the results in a list of tuple of matrix with parts
# of the answer

FGMs_index = np.zeros((nOmega,ndOmega))
for res in results:
    FGMs_index += res    

# ext = np.array([Omega_range[0]/OmegaSun, Omega_range[nOmega-1]/OmegaSun, 
#                 dlnOmegadlnr_range[0], dlnOmegadlnr_range[ndOmega-1]])

# plot
# plt.imshow(FGMs_index.transpose(), aspect="auto", interpolation="none",
#            origin="lower", extent=ext)

# plt.gca().invert_yaxis()
# plt.colorbar()
    
# plt.xlabel("$\Omega/\Omega_\odot$")
# plt.ylabel("$\partial\log\Omega / \partial \log r $")
# plt.title("$\sigma_{FGM}$ in the $\Omega$, $\partial\log\Omega / \partial r$"+
#           "space")
# plt.show()

try:
    import cpickle as pickle
except:
    import pickle

pickle.dump(FGMs_index, open("dump", "bw"))
