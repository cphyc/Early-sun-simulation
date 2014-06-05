try:
    from crunner import loop, get_param
    print("Using Cython version")
except:
    from runner import loop, get_param
    print("Using Python version")
    
import numpy as np
import sys
try:
    from multiprocessing import Pool
    multiproc = int(sys.argv[1])
    pool = Pool(processes=multiproc)
    print ("Using %d processors" % multiproc)
except:
    multiproc = False
    
from itertools import product8
	
import signal
import sys

def join_and_dump(results,FGMS_index):
    for res in results:
        FGMs_index += res    
    try:
        import cpickle as pickle
    except:
        import pickle

    pickle.dump(FGMs_index, open("dump", "bw"))

# Extension of the omega domain and the k domain
nk = 500
nOmega = 200
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
    # Create the list of parameters we want
    for nom in range(nOmega//4 + 1):
        for ndom in range(ndOmega//4 + 1):
            om_b = 4*nom
            om_e = min(om_b+4, nOmega)
            dom_b= 4*ndom
            dom_e= min(dom_b+4,ndOmega)

            params.append( [Omega_range, dlnOmegadlnr_range,
                            k_range, om_b, om_e, dom_b, dom_e] )

    # Launch nprocess threads, dump and start again
    for param_short in params[::multiproc]
        results = pool.map(loop, param_short)
        joind_and_dump(results, FGMs_index)
else:
    results = [loop([Omega_range, dlnOmegadlnr_range,
                    k_range, 0, nOmega, 0, ndOmega])]

join_and_dump(results, FGMs_index)
