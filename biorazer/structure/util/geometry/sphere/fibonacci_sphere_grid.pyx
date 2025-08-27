# distutils: language = c
# cython: boundscheck = False
# cython: wraparound = False
# cython: initializedcheck = False
# cython: cdivision = True

import numpy as np
cimport numpy as np
from libc.math cimport acos, cos, sin, pi, sqrt

def fibonacci_sphere_grid(center, double radius, int n_points, double offset=0.5):
    """
    Generate points uniformly distributed on the surface of a sphere using the Fibonacci lattice method.

    Parameters
    ----------
    center : array-like
        The (x, y, z) coordinates of the sphere's center.
    radius : float
        The radius of the sphere.
    n_points : int
        The number of points to generate on the sphere's surface.

    Returns
    -------
    points : np.ndarray(dtype=np.float64, ndim=2)
        An array of shape (n_points, 3) containing the (x, y, z) coordinates of the points on the sphere's surface.
    """
    cdef double golden_ratio = (sqrt(5.0) + 1) / 2
    cdef double x0 = center[0], y0 = center[1], z0 = center[2]
    cdef int i
    cdef double u, v, theta, phi

    cdef double[:, ::1] points = np.empty((n_points, 3), dtype=np.double)
    
    for i in range(n_points):
        u = i / golden_ratio
        u = u - <int>u  # 替代 % 1，更高效
        v = (i + offset) / (n_points - 1 + offset * 2)
        theta = 2 * pi * u
        phi = acos(1 - 2 * v)
        points[i, 0] = x0 + radius * cos(theta) * sin(phi)
        points[i, 1] = y0 + radius * sin(theta) * sin(phi)
        points[i, 2] = z0 + radius * cos(phi)
    
    return np.asarray(points)