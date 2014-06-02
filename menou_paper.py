import matplotlib.pyplot as plt
from runner import loop
import numpy as np


# Extension of the omega domain and the k domain
nOmega = 20
ndOmega = nOmega
nk = 20

# container for the results
FGMs = np.zeros((nOmega*ndOmega, 3), dtype="float64")
FGMs_index = np.zeros((nOmega,ndOmega), dtype="float64")

# run the loop
ext = loop(FGMs, FGMs_index, nOmega, ndOmega, nk)

# plot
plt.imshow(FGMs_index.transpose(), aspect="auto", interpolation="none",
           origin="lower", extent=ext)

plt.gca().invert_yaxis()
plt.colorbar()
    
plt.xlabel("$\Omega/\Omega_\odot$")
plt.ylabel("$\partial\log\Omega / \partial \log r $")
plt.title("$\sigma_{FGM}$ in the $\Omega$, $\partial\log\Omega / \partial r$"+
          "space")
plt.show()
