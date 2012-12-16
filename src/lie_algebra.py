"""
Routines for simple tasks with SE(2), such as computing the Cayley map,
its derivative, etc.

"""

from __future__ import division
import numpy as np


J = np.array([[0, 1], [-1, 0]]) # Symplectic matrix

def skew_product(v, w):
    """
    Skew product of two vectors.

    """
    return v[0]*w[1] - v[1]*w[0]


def cayley_se2(xi):
    """
    Computes the Cayley transform of an element in the Euclidian Lie
    algebra se(2).

    Parameters
    ----------

    xi : array-like (3x1)
         Lie algebra element in se(2)

    Returns
    -------

    g : array-like (3x3)
        Group element in SE(2)

    """
    omega, v, w = xi
    return np.array([
            [1-omega**2/4, -omega, v-omega*w/2],
            [omega, 1-omega**2/4, w+omega*v/2],
            [0, 0, omega**2/4+1]], dtype=np.double)/(1+omega**2/4)

def cayley_so2(omega):
    """
    Computes the Cayley transform of an element in the 2D rotation group SO(2).

    Parameters
    ----------

    xi : real number
         Angular velocity, interpreted as an element in the Lie algebra so(2)

    Returns
    -------

    g : array-like (2x2)
        Group element in SO(2)

    """
    return np.array([[1-omega**2/4, -omega],
                     [omega, 1-omega**2/4]], dtype=np.double)/(1+omega**2/4)

def delta_dcayley_inv_se2(xi, delta_xi):
    """
    Compute the directional derivative of the right-trivialized
    derivative of the Cayley map in a given direction.

    """
    omega, v, w = xi
    d_omega, d_v, d_w = delta_xi

    row1 = [omega*d_omega/2, -d_w/2+(d_omega*v+omega*d_v)/4, 
            d_v/2+(d_omega*w+omega*d_w)/4]

    return np.array([row1, [0, 0, -d_omega/2], [0, d_omega/2, 0]])
    
def dcay_inv_se2(xi):
    """
    Returns the inverse of the right-trivialized derivative of the Cayley map.

    Parameters
    ----------

    xi : array-like (3x1)
         Lie algebra element in se(2)

    Returns
    -------

    M : array-like (3x3)
        Matrix of a linear transformation on se(2)

    """
    omega, v, w = xi
    return np.array([[1+omega**2/4, 0, 0], 
                     [-w/2+omega*v/4, 1, omega/2], 
                     [v/2+omega*w/4, -omega/2, 1]]) 

def rotation(angle):
    """
    Construct a rotation matrix over a given angle.

    Parameters
    ----------

    theta : real number
            Rotation angle

    Returns
    -------

    R : array-like (2x2)
        Rotation matrix

    """
    c = np.cos(angle); s = np.sin(angle)
    return np.array([[c, -s], [s, c]])
