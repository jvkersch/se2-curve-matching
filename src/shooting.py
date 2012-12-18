"""
Solve the BVP through a shooting algorithm, using Newton iteration.

"""

import numpy as np
from lie_algebra import J
from ode_system import (integrate, coefficients_first_variation, 
                        right_boundary_condition)

import ipdb

def solve_bvp(c0, c1, m, theta_guess, full_output=False, tol=1e-8, n_iter=20):

    delta = -1; n = 0
    
    while (n < n_iter) and (abs(delta) > tol): 

        # ipdb.set_trace()

        # Integrate forward in time
        theta, omega, v, delta_theta, delta_omega, delta_v = \
            integrate(c0, c1, theta_guess, m)
        
        # Calculate adjustment for angle as the fraction numer/denom
        numer = (m*(1+omega[-1]**2/4)*omega[-1] + 
                 omega[-1]/4*np.dot(v[-1, :], v[-1, :]) + 
                 np.dot(v[-1, :], np.dot(J, c0[-1, :])) - 
                 omega[-1]/2*np.dot(v[-1, :], c0[-1, :]))
            
        c, d = coefficients_first_variation(omega[-1], v[-1, :], c0[-1, :], m)
        denom = c*delta_omega[-1] + np.dot(d, delta_v[-1, :])
         
        # Update initial guess
        n += 1
        delta = numer/denom
        theta_guess -= delta        
    
    out = theta, omega, v, delta_theta, delta_omega, delta_v
    if full_output:
        res = right_boundary_condition(theta[-1], omega[-1], v[-1, :], 
                                       c0[-1, :], c1[-1, :], m)
        out += (res, n)

    return out


