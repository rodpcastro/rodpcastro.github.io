import numpy as np
from twodubem.geometry import Polygon
from twodubem._internal import tozero

class FloatingCylinder(Polygon):
    """Floating cylinder.

    Parameters
    ----------
    radius : float
        Cylinder radius.
    number_of_elements : int
        Number of elements along the circumference.
    """

    def __init__(self, radius, number_of_elements):
        self.radius = radius
        self.number_of_elements = number_of_elements
        self._set_vertices()
        self._set_elements()
        self._set_sides_properties()

    def _set_vertices(self):
        t = np.linspace(0.0, -np.pi, self.number_of_elements+1)
        x = self.radius * np.cos(t)
        y = self.radius * np.sin(t)
        vertices = np.column_stack((x, y))
        self.vertices = tozero(vertices)

    def show(self, filename=''):
        """Display a graphical representation of the boundary."""

        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon as PolygonPatch
        
        fig, ax = plt.subplots()
        
        # Parameters for framing the mesh.
        xmax = 1.2*self.radius
        ymax = 0.2*self.radius
        xmin = ymin = -xmax
        
        ax.spines[:].set_visible(False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.margins(0)
        
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        
        ax.set_aspect('equal')
        
        ax.set_facecolor('skyblue')
        
        # Circular section.
        cpatch_vertices = np.zeros((self.vertices.shape[0]+1,
                                    self.vertices.shape[1]))
        cpatch_vertices[:-1] = self.vertices
        cpatch_vertices[-1] = self.vertices[0]
        cpatch = PolygonPatch(cpatch_vertices, color='white')
        
        # Region above free-surface.
        rpatch_vertices = np.array([[xmin, 0.0 ],
                                    [xmax, 0.0 ],
                                    [xmax, ymax],
                                    [xmin, ymax]])
        rpatch = PolygonPatch(rpatch_vertices, color='white')
        
        ax.add_patch(cpatch)
        ax.add_patch(rpatch)
        ax.plot(
            self.vertices[:, 0],
            self.vertices[:, 1],
            'r-',
            marker='o',
            markersize=3,
            markerfacecolor='b',
            markeredgecolor='k',
            linewidth=2,
        )

        if filename:
            plt.savefig(filename, bbox_inches='tight', pad_inches=0)

        plt.show()