"""
Integrate curve matching equations and compute boundary conditions.

Do not use these routines directly. Instead, use ... in 

"""
import math
import numpy as np
import scipy.optimize as so

from lie_algebra import rotation, cayley_so2, J, skew_product


def projected_momentum(omega, v, pi, p, c_pt):
    """
    Compute the projection of an se(2)-momentum, acted upon by the 
    adjoint of the d-Cayley-inverse map, using the optimized representation
    given in the paper.

    Parameters
    ----------

    omega, v : scalar, ndarray
               Angular and linear velocity.

    pi, p : scalar, ndarray
            Angular and linear momentum.

    c_pt : ndarray
           Basepoint in R2.

    Returns
    -------

    proj : scalar
           Projection onto stabilizer of ``c_pt``.

    """

    return ((1 + omega**2/4)*pi - 1./2*skew_product(p, v) + 
            omega/4*np.dot(p, v) + skew_product(p, c_pt) - 
            omega/2*np.dot(p, c_pt))


def compute_linear_displacement(c0_current, c1_current, theta):
    """
    Solve constraint equation for linear displacement.

    Parameters
    ----------

    c0_current, c1_current : ndarray
                             Corresponding points on matching curves.
                
    theta : scalar
            Angle matching ``c0_current`` to ``c0_current``.

    """
    return c1_current - np.dot(rotation(theta), c0_current)


def compute_linear_velocity(R1, omega1, c0_points, c1_points):
    """
    Solve time-derivative of constraint equation for linear velocity.

    Parameters
    ----------

    R1 : ndarray
         Rotation matrix matching the points in step ``k`` of the algorithm.

    omega1 : scalar
             Angular velocity to update rotation from step ``k`` to ``k+1``.

    c0_points, c1_points: array_like
                          Tuples of ndarrays, representing the points to be
                          matched at steps ``k`` and ``k+1``.

    Returns
    -------

    v : ndarray
        Linear velocity satisfying the (time-derivative of the) constraint.

    """
    
    c0_current, c0_next = c0_points
    c1_current, c1_next = c1_points

    b = (np.dot(R1.T, c1_next - c1_current) + c0_current - 
         np.dot(cayley_so2(omega1), c0_next))
    B = np.array([[1, omega1/2], [-omega1/2, 1]])
    
    return np.dot(B, b)


def integrate_one_step(R1, omega0, v0, c0_points, c1_points, m):
    """
    Perform one step in the integration algorithm for the curve matching 
    equations.

    Parameters
    ----------

    R1 : ndarray
         Rotation matrix matching the points in step ``k`` of the algorithm.

    omega0 : scalar
             Angular velocity from previous iteration of algorithm.

    v0 : ndarray
         Linear velocity from previous iteration of algorithm.

    c0_points, c1_points: array_like
                          Tuples of ndarrays, representing the points to be
                          matched at steps ``k`` and ``k+1``.
                          
    m : scalar
        Relative weight for angular deformation.

    Returns
    -------

    R2 : ndarray
         Rotation matrix matching the points in step ``k+1`` of the algorithm.

    omega1, v1 : scalar, ndarray
                 Components of angular and linear velocity in step ``k+1``.

    """

    # Unpack points
    # 
    # c0_current has already been matched to c1_current,
    # c0_next is the one we're trying to match to c1_next
    c0_current, c0_next = c0_points
    c1_current, c1_next = c1_points

    rhs = projected_momentum(-omega0, -v0, m*omega0, v0, c0_current)
    def optimization_function(omega1):
        # Scalar function to find next omega
        v1 = compute_linear_velocity(R1, omega1, c0_points, c1_points)
        return projected_momentum(omega1, v1, m*omega1, v1, c0_current) - rhs

    # Determine components of Lie algebra update element
    omega1 = so.newton(optimization_function, omega0)
    v1 = compute_linear_velocity(R1, omega1, c0_points, c1_points)

    # Determine next group matching element
    # Only the rotational part is needed, since the translational part 
    # can be determined from the constraints, if necessary.
    R1 = np.dot(R1, cayley_so2(omega1))

    return R1, omega1, v1


def coefficients_first_variation(omega, v, c_pt, m):
    """
    Helper function to compute coefficients of first-variation equation.

    """
    c = m*(1+3./4*omega**2) + np.dot(v, v)/4 - np.dot(v, c_pt)/2
    d = omega/2*v + np.dot(J, c_pt) - omega/2*c_pt
    return c, d


def compute_first_variation(omega0, v0, omega1, v1, R1, 
                            delta_omega0, delta_v0, delta_theta1,
                            c0_points, c1_points, m):
    """
    Solve the first-variation equations. This amounts to solving a 
    (complicated) system of linear equations.

    Parameters
    ----------

    omega0, v0 : scalar, ndarray
                 Angular and linear velocity at step ``k-1``.

    omega1, v1 : scalar, ndarray
                 Angular and linear velocity at step ``k``.

    R1 : ndarray
         Rotation matrix at step ``k``.

    delta_omega0, delta_v0: scalar, ndarray
                            Components of first variation at step ``k-1``.

    delta_theta1 : scalar
                   First variation of angle at step ``k``.

    c0_points, c1_points: array_like
                          Tuples of ndarrays, representing the points to be
                          matched at steps ``k`` and ``k+1``.
                          
    m : scalar
        Relative weight for angular deformation.

    Returns
    -------

    delta_omega1, delta_v1, delta_theta2 : scalar, array_like, scalar

    """

    c0_current, c0_next = c0_points
    c1_current, c1_next = c1_points

    c_plus, d_plus = coefficients_first_var(omega1, v1, c0_current, m)
    c_neg, d_neg = coefficients_first_var(-omega0, -v0, c0_current, m)

    b = (np.dot(R1.T, c1_next - c1_current) + c0_current - 
         np.dot(cayley_so2(omega1), c0_next))

    B = np.array([[1, omega1/2], [-omega1/2, 1]])
    A = np.linalg.inv(B)

    C = np.dot(J, b/2 + np.dot(A, c0_next))
    D = np.dot(B, np.dot(R1.T, np.dot(J, c1_next - c1_current)))

    delta_omega1 = 1./(c_plus + np.dot(d_plus, C))*(
        c_neg*delta_omega0 + np.dot(d_neg, delta_v0) - 
        np.dot(d_plus, D)*delta_theta1)

    delta_v1 = C*delta_omega1 + D*delta_theta1

    delta_theta2 = (delta_theta1 + 
                    (1-omega1**2/4)/(1+omega1**2/4)*delta_omega1)
    
    return delta_omega1, delta_v1, delta_theta2


def integrate(c0, c1, theta0, m):
    """
    Integrate curve matching equations for two given equations.

    Parameters
    ----------

    c0, c1 : ndarray
             Curves to match. Must have the same number of points

    theta0 : scalar
             Initial value for matching angle.

    m : scalar
        Relative weight for angular deformation.
       
    Returns
    -------

    theta : ndarray
            Array of angles along the curve.

    omega, v : ndarray, ndarray
               Components of angular and linear velocity along the curve.

    delta_theta, delta_omega, delta_v : ndarray, ndarray, ndarray
                                        Components of first variation along 
                                        the curve.

    """

    # Book-keeping
    numpoints = c0.length()

    theta = np.zeros(numpoints)
    omega = np.zeros(numpoints-1)
    v = np.zeros((numpoints-1, 2))

    delta_theta = np.zeros(numpoints)
    delta_omega = np.zeros(numpoints-1)
    delta_v = np.zeros((numpoints-1, 2))

    # Set up initial conditions
    R0 = rotation(theta0); R1 = R0
    omega0 = 0
    v0 = compute_linear_velocity(R0, omega0, 
                                 (c0[0, :], c0[1, :]), 
                                 (c1[0, :], c1[1, :]))

    delta_theta0 = 1
    delta_theta1 = 1
    
    delta_omega0 = 0
    delta_v0 = np.dot(R0.T, np.dot(J, c1[1]-c1[0]))*delta_theta0

    # Record initial values
    theta[0] = theta0; theta[1] = theta0
    omega[0] = omega0 # Zero
    v[0, :] = v0

    delta_theta[0] = delta_theta0
    delta_theta[1] = delta_theta1
    delta_omega[0] = delta_omega0
    delta_v[0, :] = delta_v0

    # Main integration loop 
    for s in xrange(1, numpoints-1):

        c0_points = (c0[s, :], c0[s+1, :])
        c1_points = (c1[s, :], c1[s+1, :])

        # Compute next element of solution trajectory
        R2, omega1, v1 = integrate_one_step(R1, omega0, v0, 
                                            c0_points, c1_points, m)

        # Compute next step of first-variation equation
        delta_omega1, delta_v1, delta_theta2 = \
            compute_first_variation(omega0, v0, omega1, v1, R1, 
                                    delta_omega0, delta_v0, delta_theta1,
                                    c0_points, c1_points, m)

        # Record computed values
        omega[s] = omega1
        theta[s+1] = math.atan2(R1[1,0], R1[0,0]) 

        v[s, :] = v1

        delta_omega[s] = delta_omega1
        delta_theta[s+1] = delta_theta2 

        delta_v[s, :] = delta_v1

        # Move along
        R1 = R2
        omega0 = omega1
        v0 = v1

        delta_theta1 = delta_theta2
        delta_omega0 = delta_omega1
        delta_v0 = delta_v1

    return theta, omega, v, delta_theta, delta_omega, delta_v


def right_boundary_condition(theta, omega, v, c0, c1, m):
    """
    Residual of boundary condition at right terminal end of curve.

    Parameters
    ----------

    theta, omega, v : scalar, scalar, ndarray
                      Terminal matching angle, angular and linear velocity.

    c0_end, c1_end : ndarray, ndarray
                     Terminal points of the curves.

    m : scalar
        Relative weight for angular deformation.

    Returns
    -------

    res : scalar
          Residual for right terminal boundary condition.

    """
    return ( m*(1+omega**2/4)*omega + omega/4*np.dot(v, v) + 
             np.dot(v, np.dot(J, c0)) - omega/2*np.dot(v, c0) )
