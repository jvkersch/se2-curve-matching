import numpy as np

from scipy.linalg import expm, logm

from lie_algebra import rotation
from curves import DiscreteCurve
from matching import GroupTrajectory

class LiegroupInterpolationMethod:

    def __init__(self, c0, c1, g):

        self.N = c0.length()
        self.c0_3_raw = np.empty((self.N, 3))
        self.c1_3_raw = np.empty((self.N, 3))
        self.xi = np.empty((3, 3, self.N))

        # Copy c0, c1 into Nx3 arrays with ones in the last column
        self.c0_3_raw[:, 0:2] = c0.points
        self.c1_3_raw[:, 0:2] = c1.points
        self.c0_3_raw[:, 2] = np.ones(self.N)
        self.c1_3_raw[:, 2] = np.ones(self.N)

        # Compute infinitesimal transformations
        for k in xrange(0, self.N):
            gmat = np.empty((3, 3))
            gmat[0:2, 0:2] = rotation(g.theta[k])
            gmat[0:2, 2] = g.x[k, :]
            gmat[2, 2] = 1

            # Obtain Lie algebra element
            self.xi[:, :, k] = logm(gmat)

    def transform(self, epsilon):

        c0_transformed = DiscreteCurve(N=self.N)
        for k in xrange(0, self.N):
            g = expm(epsilon*self.xi[:, :, k])
            c0_transformed[k, :] = np.dot(g, self.c0_3_raw[k, :])[0:2]

        return c0_transformed
