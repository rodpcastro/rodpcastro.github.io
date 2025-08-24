import numpy as np
from twodubem.geometry import Polygon
from twodubem._internal import tozero


class Body(Polygon):
    """Floating body."""

    def __init__(self):
        self.dofs = {}
        self.M = np.empty(0)
        self.C = np.empty(0)

    def _set_masks(self):
        self._mask_body = np.zeros(self.number_of_elements, dtype=np.bool)
        self._mask_lid = np.zeros(self.number_of_elements, dtype=np.bool)
        self._mask_body[:self.number_of_body_elements] = True
        self._mask_lid[self.number_of_body_elements:] = True

    def add_degree_of_freedom(self, name, cg=np.array([0.0, 0.0])):
        name_ = name.strip().lower()
        midpoints = self.midpoints[self._mask_body]
        normals = self.normals[self._mask_body]

        if name_ == 'sway':
            dof = normals[:, 0]
        elif name_ == 'heave':
            dof = normals[:, 1]
        elif name_ == 'roll':
            x = midpoints - cg
            dof = x[:, 0] * normals[:, 1] - x[:, 1] * normals[:, 0]
        else:
            raise ValueError('Invalid degree of freedom')

        self.dofs[name_] = dof
    
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


class Cylinder(Body):
    """Floating cylinder.

    Parameters
    ----------
    radius : float
        Cylinder radius.
    number_of_elements : int
        Number of elements along the cylinder's circumference.
    """

    def __init__(self, radius, number_of_elements):
        super().__init__()
        self.radius = radius
        self.number_of_body_elements = number_of_elements
        self._set_vertices()
        self._set_elements()
        self._set_masks()
        self._set_sides_properties()

    def _set_vertices(self):
        # Cylinder vertices.
        t = np.linspace(0.0, -np.pi, self.number_of_body_elements+1)
        x = self.radius * np.cos(t)
        y = self.radius * np.sin(t)
        cylinder_vertices = tozero(np.column_stack((x, y)))

        # Lid vertices.
        elength = self.radius * np.pi / self.number_of_body_elements
        nlid = np.floor(2 * self.radius / elength).astype(int)
        self.number_of_lid_elements = nlid
        lid_vertices = np.linspace(cylinder_vertices[-1],
                                   cylinder_vertices[0],
                                   self.number_of_lid_elements+1)

        vertices = np.vstack((cylinder_vertices, lid_vertices[1:]))
        self.vertices = vertices

    def set_inertia_matrix(self, rho=1.0):
        m = rho * 0.5 * np.pi * self.radius**2
        ndofs = len(self.dofs)
        self.M = np.zeros((ndofs, ndofs))
        for i, dof in enumerate(self.dofs.keys()):
            if dof in ['sway', 'heave']:
                self.M[i, i] = m
            elif dof == 'roll':
                self.M[i, i] = 0.5 * m * self.radius**2

    def set_stiffness_matrix(self, rho=1.0, g=9.81):
        ndofs = len(self.dofs)
        self.C = np.zeros((ndofs, ndofs))
        for i, dof in enumerate(self.dofs.keys()):
            if dof == 'heave':
                self.C[i, i] = 2 * self.radius * rho * g
            elif dof in ['roll', 'sway']:
                self.C[i, i] = 0.0