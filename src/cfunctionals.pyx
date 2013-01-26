# cython: profile=False

import numpy as np
cimport numpy as np
cimport cython 

DTYPE=np.double
ctypedef np.double_t DTYPE_t

def energy_functional_slow(thetas, c0, c1, m, h):
    """
    Very slow version of the energy functional. To be used for comparison 
    purposes.

    """
    import math
    from lie_algebra import cayley_so2, rotation

    E = 0 
    bound = thetas.shape[0] - 1

    for k in xrange(0, bound):
        omega = 2*math.tan((thetas[k+1] - thetas[k])/2)
        
        R = rotation(thetas[k])

        b = (np.dot(R.T, c1[k+1,:] - c1[k,:]) + c0[k,:] - 
             np.dot(cayley_so2(omega), c0[k+1,:]))

        E += (m*omega**2 + (1+omega**2/4)*np.dot(b, b))/(2*h)
            
    return E


@cython.cdivision(True)
@cython.boundscheck(False) 
@cython.wraparound(False)
cpdef double energy_functional(np.ndarray[DTYPE_t, ndim=2] X, 
                               np.ndarray[DTYPE_t, ndim=2] c0, 
                               np.ndarray[DTYPE_t, ndim=2] c1, 
                               double m, double h):
    """
    Optimized version of the energy functional.

    """
    
    cdef int k, bound

    cdef double E = 0
    cdef double omega

    # Auxiliary coordinates
    cdef double vx, vy, c0x, c0y, x0, y0, x1, y1, bx, by
    cdef double C, S, factor, denominator

    bound = c0.shape[0] - 1

    for k from 0 <= k < bound:
        # Assign array elements to coordinates individually,
        # since Cython is not very good at handling slicing.
        x0 = X[k, 0]
        y0 = X[k, 1]

        x1 = X[k+1, 0]
        y1 = X[k+1, 1]

        # Components of the cos and sin of the difference of angles 
        C = x1*x0 + y1*y0
        S = y1*x0 - x1*y0
        omega = 2*S/(1+C)

        vx = c1[k+1, 0] - c1[k, 0]
        vy = c1[k+1, 1] - c1[k, 1]
        c0x = c0[k+1, 0]
        c0y = c0[k+1, 1]

        factor = 1 - omega*omega/4
        denominator = 1 + omega*omega/4

        # b quantity proportional to the velocity
        bx = c0[k, 0] - (factor*c0x - omega*c0y)/denominator
        by = c0[k, 1] - (omega*c0x + factor*c0y)/denominator
        
        bx +=  x0*vx + y0*vy
        by += -y0*vx + x0*vx

        E += (m*omega*omega + (1+omega*omega/4)*(bx*bx + by*by))/(2*h)
            
    return E
