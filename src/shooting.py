"""
Solve the BVP through a shooting algorithm, using Newton iteration.

"""

import numpy as np
from matching import (integrate, energy, J, 
                      coefficients_first_var, right_boundary_condition)


class ShootingMethod:
    def __init__(self, tol=1e-8, n_iter=20):

        self.tol = tol
        self.n_iter = n_iter

    def iterate(self, c0, c1, m, theta_guess, full_output=False):

        delta = -1; n = 0
    
        while (n < self.n_iter) and (abs(delta) > self.tol): 
            # Integrate forward in time
            theta, omega, v, delta_theta, delta_omega, delta_v = \
                integrate(c0, c1, theta_guess, m, full_output=True)
        
            # Calculate adjustment for angle as the fraction numer/denom
            numer = (m*(1+omega[-1]**2/4)*omega[-1] + 
                     omega[-1]/4*np.dot(v[-1, :], v[-1, :]) + 
                     np.dot(v[-1, :], np.dot(J, c0[-1, :])) - 
                     omega[-1]/2*np.dot(v[-1, :], c0[-1, :]))
            
            c, d = coefficients_first_var(omega[-1], v[-1, :], c0[-1, :], m)
            denom = c*delta_omega[-1] + np.dot(d, delta_v[-1, :])
         
            # Update initial guess
            n += 1
            delta = numer/denom
            theta_guess -= delta        
    
        out = theta, omega, v, delta_theta, delta_omega, delta_v
        if full_output:
            res = right_boundary_condition(theta, omega, v, c0, c1, m)
            out += (res, n)

        return out


def sweep_theta_range(c0, c1, h, m, N=100):
    """
    Sweep across a range of initial guesses for theta and record energy of 
    corresponding minimizer, obtained through shooting.

    """
    energies = np.empty((N, 4))
    s = ShootingMethod()
    
    for n in xrange(0, N):
        theta_guess = 2*np.pi*n/N
        print "%d: theta_guess = %f, " % (n, theta_guess),

        _, omega, v, _, _, _, res, n_iter = \
            s.iterate(c0, c1, m, theta_guess, full_output=True)
        E = energy(omega, v, h, m)
        energies[n, :] = theta_guess, E, res, n_iter

        print "E = %e, res = %e, n_iter = %d" % (E, res, n_iter)

    return energies
