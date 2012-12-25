import numpy as np

import math
from lie_algebra import cayley_so2, rotation




def energy_functional(thetas, *args):
    
    c0, c1, m, h = args
    
    E = 0
    
    #theta_current = theta0
    for k in xrange(0,thetas.shape[0]-1):
        #theta_next = thetas[k]
        omega = 2*math.tan((thetas[k+1] - thetas[k])/2)
        
        R = rotation(thetas[k])
        #c0_points = (c0[k, :], c0[k+1, :])
        #c1_points = (c1[k, :], c1[k+1, :])
        
        # v = compute_linear_velocity(R, omega, c0_points, c1_points)
        b = (np.dot(R1.T, c1[k+1,:] - c1[k,:]) + c0[k,:] - 
             np.dot(cayley_so2(omega1), c0[k+1,:]))


        E += (m*omega**2 + (1+omega**2/4)*np.dot(b, b))/(2*h)
        
        #theta_current = theta_next
        
        # print omega, v
    
    return E
