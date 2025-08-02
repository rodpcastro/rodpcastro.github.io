import numpy as np
from twodubem.geometry import Polygon


class FloatingCylinder(Polygon):
    """Fluid domain of a floating cylinder.

    Domain size is based on lmax.
    Mesh size is based on lmin.

    Parameters
    ----------
    lmin : float
        Minimum wave length.
    lmax : float
        Maximum wave length.

    The diagram below describes the number of elements on each boundary.

    ×--nf--×    ×--nf--×
    |      |_nc_|      |
    nd                 nd
    |                  |
    ×--------nb--------×

    nf: number of free surface elements on each side of the floating cylinder.
    nd: number of depth elements on each side of the floating cylinder.
    nb: number of bottom elements.
    nc: number of elements on the floating cylinder.
    """

    def __init__(self, lmin, lmax):
        self.lmin = lmin
        self.lmax = lmax

        self.radius = 1.0
        self.water_depth = np.max([2.0*self.radius, self.radius + lmax])
        self.radiation_radius = np.max([2.0*self.radius, self.radius + 2*lmax])

        self.msize = self.lmin / 16  # Reference element size.
        self.fr = 1.02  # Reference ratio for free surface geometric progression.
        self.hr = 1.05  # Reference ratio for depth geometric progression.

        self._set_mesh_parameters()
        self._set_elements()
        self._set_sides_properties()
        self._set_boundary_masks()
        self._set_boundary_orientation()
        self._set_boundary_size()
    
    def _set_mesh_parameters(self):
        c0 = np.array([0.0, 0.0])  # Cylinder's center
        xv = np.array([1.0, 0.0])  # x unit vector
        yv = np.array([0.0, 1.0])  # y unit vector
        R = self.radius
        h = self.water_depth
        s = self.radiation_radius

        # Domain corners.
        v1 = c0 - R * xv
        v2 = v1 - (s-R) * xv
        v3 = v2 - h * yv
        v4 = v3 + 2*s * xv
        v5 = v4 + h * yv
        v6 = v5 - (s-R) * xv

        # Free surface vertices.
        xv, self.nf, self.fr = geometric_progression(R, s, self.msize, self.fr)

        free_surface_n = np.zeros((self.nf, 2))
        free_surface_n[0] = v1
        free_surface_n[1:, 0] = -xv[1:-1]

        free_surface_p = np.zeros((self.nf, 2))
        free_surface_p[0] = v5
        free_surface_p[1:, 0] = xv[-2:0:-1]
        
        # Depth vertices.
        depth_mesh_size = self.msize * self.fr**self.nf
        yv, self.nd, self.hr = geometric_progression(0.0, h, depth_mesh_size, self.hr)

        depth_n = np.zeros((self.nd, 2))
        depth_n[:, 0] = v2[0]
        depth_n[:, 1] = -yv[:-1]

        depth_p = np.zeros((self.nd, 2))
        depth_p[:] = v4
        depth_p[1:, 1] = -yv[-2:0:-1]

        # Bottom vertices.
        bottom_mesh_size = depth_mesh_size * self.hr**self.nd
        self.nb = np.ceil(2*s/bottom_mesh_size).astype(int)
        bottom = np.linspace(v3, v4, self.nb, endpoint=False)

        # Cylinder vertices.
        self.nc = np.max([6, np.pi*self.radius/self.msize]).astype(int)
        theta = np.linspace(0.0, -np.pi, self.nc, endpoint=False)
        cylinder = R * np.column_stack((np.cos(theta), np.sin(theta)))

        # Combine all vertices, in order.
        self.vertices = np.concatenate(
            (
                free_surface_n,
                depth_n,
                bottom,
                depth_p,
                free_surface_p,
                cylinder,
                free_surface_n[:1],
            )
        )

    def _set_boundary_masks(self):
        self.free_surface = np.zeros(self.number_of_elements, dtype=np.bool)
        self.depth = np.zeros(self.number_of_elements, dtype=np.bool)
        self.bottom = np.zeros(self.number_of_elements, dtype=np.bool)
        self.cylinder = np.zeros(self.number_of_elements, dtype=np.bool)

        n1 = 0
        n2 = n1 + self.nf
        n3 = n2 + self.nd
        n4 = n3 + self.nb
        n5 = n4 + self.nd
        n6 = n5 + self.nf

        self.free_surface[n1:n2] = True
        self.free_surface[n5:n6] = True
        self.depth[n2:n3] = True
        self.depth[n4:n5] = True
        self.bottom[n3:n4] = True
        self.cylinder[n6:] = True


def geometric_progression(x1, x2, l0, rt):
    s = x2 - x1

    n = np.ceil(np.log(s*(rt-1)/l0 + 1) / np.log(rt)).astype(int)

    # Adjust geometric progression ratio to match calculated n.
    f = lambda r: r**n + s * (1-r) / l0 - 1
    fp = lambda r: n * r**(n-1) - s / l0

    rtn = newton_raphson(f, fp, rt)

    # Generate points.
    xv = np.empty(n+1)
    xv[0] = x1

    for i in range(n):
        xv[i+1] = xv[i] + l0 * rtn**i

    return xv, n, rtn


def newton_raphson(f, fp, x0):
    for i in range(20):
        x1 = x0 - f(x0) / fp(x0)
        x0 = x1
        if f(x1) < 1.0e-8:
            break

    return x1