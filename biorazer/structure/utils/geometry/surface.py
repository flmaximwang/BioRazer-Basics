import numpy as np
from scipy.spatial import KDTree
from .sphere import fibonacci_sphere_grid


def fibonacci_surface_grid(
    coords, vdw_radii, probe_radius=1.4, n_points=1000, offset=0.5, shrink=0.01
):
    """
    Generate points uniformly distributed on the surface of spheres centered at given coordinates.

    Parameters
    ----------
    coords : array-like, shape (n, 3)
        The (x, y, z) coordinates of the sphere centers.
    vdw_radii : array-like, shape (n,)
        The van der Waals radii corresponding to each coordinate.
    probe_radius : float, optional
        The radius of the probe sphere, by default 1.4.
    n_points : int, optional
        The number of points to generate on each sphere's surface, by default 1000.
    offset : float, optional
        An offset to avoid placing points directly at the poles, by default 0.5.
    shrink : float, optional
        A small value to shrink the effective radius for point generation, by default 0.01.
        If no shrink is applied, some points may be removed undesirably because of numerical inaccuracies.

    Returns
    -------
    points : np.ndarray
        An array of shape (n * n_points, 3) containing the (x, y, z) coordinates of the points on the spheres' surfaces.
    """
    selector = vdw_radii > 0
    real_coords = coords[selector]
    real_vdw_radii = vdw_radii[selector]
    surface_points = np.zeros((len(real_coords) * n_points, 3))
    counter = 0
    for coord, vdw_radius in zip(real_coords, real_vdw_radii):
        surface_points[counter : counter + n_points] = fibonacci_sphere_grid(
            center=coord,
            radius=vdw_radius + probe_radius,
            n_points=n_points,
            offset=offset,
        )
        counter += n_points

    # Remove points that are too close to other atoms
    delete_mask = np.ones(surface_points.shape[0], dtype=bool)
    surface_kdtree = KDTree(surface_points)
    delete_indices = set()
    for coord, vdw_radius in zip(real_coords, real_vdw_radii):
        delete_indices.update(
            surface_kdtree.query_ball_point(coord, r=probe_radius + vdw_radius - shrink)
        )
    delete_mask[list(delete_indices)] = False
    surface_points = surface_points[delete_mask]
    return surface_points
