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
    nprocess = int(sys.argv[1])
    pool = Pool(processes=nprocess)
    print ("Using %d processors" % nprocess)
except:
    nprocess = 0
    
from itertools import product
	
import signal
import sys

def join_and_dump(results, FGMs_index, last, nprocess):
    for res in results:
        FGMs_index += res    
    try:
        import cpickle as pickle
    except:
        import pickle
    print("Dumping at %d with %d processor%s." % (last, nprocess, "s"*(nprocess>1)))
    pickle.dump((last, nprocess, FGMs_index), open("dump", "bw"))

 
def pop_n(_in, n, first=0):
    while first <= len(_in):
        yield first, _in[first:][:n]
        first += n
    return StopIteration

        
# Extension of the omega domain and the k domain
nk = 5
nOmega = 25
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


FGMs_index = np.zeros((nOmega, ndOmega))

print ("Simulation with %d processors" % nprocess)
print ("k_range : {}\nOmega_range : {}\ndlnOmegadlnr_range : {}".format(k_range, Omega_range, dlnOmegadlnr_range))
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
print("Looping'n'dumping ...")
for last, next_params in pop_n(params, max(nprocess,1)):
    if nprocess > 0:
        results = pool.map(loop, next_params)
    else:
        results = [loop(p) for p in next_params]
        
    join_and_dump(results, FGMs_index, last, max(nprocess,1))
print("Looped !")
