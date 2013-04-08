import numpy as np

from ode_system import integrate, compute_linear_displacement
from shooting import solve_bvp


class MatchingProblem:

    def __init__(self, c0, c1, m):

        if c0.length() != c1.length():
            raise ValueError('Curves must have the same number of points.')

        self.num_points = c0.length()

        self.c0 = c0
        self.c1 = c1

        self.h = c0.h
        self.m = m

    def integrate(self, theta0):
        """
        Integrate curve matching equations along the curve. Note that this 
        trajectory does not necessarily satisfy the natural boundary 
        conditions. For this, use `match` instead.

        Parameters
        ----------

        theta0 : scalar
                 Initial value for angle

        Returns
        -------

        g : GroupTrajectory
            Trajectory in SE(2) mapping c0 to c1. 

        """
        theta, omega, v, delta_theta, delta_omega, delta_v = \
            integrate(self.c0, self.c1, theta0, self.m)

        g = GroupTrajectory(self.num_points, self.h, self.m)
        g.theta = theta 
        g.x = self.compute_displacements(theta)

        g.omega = omega
        g.v = v

        g.delta_omega = delta_omega
        g.delta_v = delta_v
        
        return g

    def match(self, theta_guess):
        """
        Matching...

        """
        theta, omega, v, delta_theta, delta_omega, delta_v, res, n_iter = \
            solve_bvp(self.c0, self.c1, self.m, theta_guess, full_output=True)

        g = GroupTrajectory(self.num_points, self.h, self.m)
        g.theta = theta 
        g.x = self.compute_displacements(theta)

        g.omega = omega
        g.v = v

        g.delta_omega = delta_omega
        g.delta_v = delta_v
        
        return g, res, n_iter

    def sweep_theta_range(self, N=100, verbose=True):
        """
        Sweep across a range of initial guesses for theta and record energy of 
        corresponding minimizer, obtained through shooting.
        
        """
        energies = np.empty((N, 4))
    
        for n in xrange(0, N):
            theta_guess = 2*np.pi*n/N

            g, res, n_iter = self.match(theta_guess)
            E = g.energy()

            energies[n, :] = theta_guess, E, res, n_iter

            if verbose: 
                print "%d: theta_guess = %f, theta_final = %f, " \
                    "E = %e, res = %e, n_iter = %d" % \
                    (n, theta_guess, g.theta[0], E, res, n_iter)

        return energies

    def compute_displacements(self, theta):
        """
        Helper function...

        """
        # Fill in linear displacements
        x = np.empty((self.num_points, 2))
        for k in xrange(0, self.num_points):
            x[k, :] = compute_linear_displacement(self.c0[k, :], 
                                                  self.c1[k, :], 
                                                  theta[k])

        return x


class GroupTrajectory:

    def __init__(self, N, h, m):

        # Placeholders for curve variables
        self.theta = None
        self.x     = None

        self.omega = None
        self.v     = None
        
        self.delta_omega = None
        self.delta_v     = None

        self.N = N

        self.h = h
        self.m = m

    def energy(self):
        """
        Energy of this curve.
        
        """
        return ( self.m*np.sum(self.omega*self.omega) + 
                 np.sum(self.v*self.v) )/(2*self.h)
