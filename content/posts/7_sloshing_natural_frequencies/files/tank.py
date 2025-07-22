import numpy as np
from twodubem.geometry import Polygon


class RectangularTank(Polygon):
    """Rectangular tank.

    Parameters
    ----------
    width : float
        Tank's width.
    depth : float
        Tank's depth.
    number_of_width_elements : int
        Number of elements on bottom of the tank.
    number_of_depth_elements : int
        Number of elements on the vertical sides of the tank.
    number_of_free_surface_elements : int
        Number of elements on the tank free surface.
    """

    def __init__(
        self,
        width,
        depth,
        number_of_width_elements,
        number_of_depth_elements,
        number_of_free_surface_elements,
    ):
        self.width = width
        self.depth = depth
        self.number_of_width_elements = number_of_width_elements
        self.number_of_depth_elements = number_of_depth_elements
        self.number_of_free_surface_elements = number_of_free_surface_elements
        self._set_vertices()
        self._set_elements()
        self._set_sides_properties()
        self._set_boundary_orientation()
        self._set_boundary_size()

    def _set_vertices(self):
        w = self.width
        h = self.depth
        nx = self.number_of_width_elements
        ny = self.number_of_depth_elements
        nz = self.number_of_free_surface_elements

        p0 = np.array([-w/2.0, 0.0])
        n = 2*ny + nx + nz
        self.vertices = np.empty((n + 1, 2))

        for side in range(4):
            if side == 0:
                # Left side.
                p1 = p0 + np.array([0.0, -h])
                i1 = ny
                self.vertices[0:i1] = np.linspace(p0, p1, ny, endpoint=False)
            elif side == 1:
                # Bottom side.
                p2 = p1 + np.array([w, 0.0])
                i2 = ny + nx
                self.vertices[i1:i2] = np.linspace(p1, p2, nx, endpoint=False)
            elif side == 2:
                # Rigth side.
                p3 = p2 + np.array([0.0, h])
                i3 = ny + nx + ny
                self.vertices[i2:i3] = np.linspace(p2, p3, ny, endpoint=False)
            elif side == 3:
                # Top side (free surface).
                p4 = p3 + np.array([-w, 0.0])
                i4 = ny + nx + ny + nz
                self.vertices[i3:i4] = np.linspace(p3, p4, nz, endpoint=False)

        # The last vertex must be equal to the first to form a closed boundary.
        self.vertices[-1] = self.vertices[0]