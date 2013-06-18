"""
Utility functions.

"""

import numpy as np

def separate_energies(energies, tol=1e-5):
    """
    Separate energies according to whether convergence was attained.

    """

    E_convergent = []
    E_divergent  = []

    for energy in energies:
        theta, e, res, n_iter = energy
        if abs(res) > 1e-5:
            E_divergent.append(energy)
        else:
            E_convergent.append(energy)
        
    E_convergent = np.array(E_convergent)
    E_divergent  = np.array(E_divergent)

    return E_convergent, E_divergent
