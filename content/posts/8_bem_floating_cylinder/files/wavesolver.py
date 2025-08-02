import numpy as np
from twodubem.solver import Solver
from twodubem.laplace import Laplace


class RadiationSolver(Solver):

    def __init__(self, boundaries):
        self.boundary = boundaries
        self.green = Laplace()
        self.method = 'constant'
        self.g = 9.81  # Acceleration of gravity
        self.rho = 1.0  # Water density

    def solve(self, w):
        """Solve the radiation potential of oscilation in heave."""

        k = w**2/self.g  # For infinite depth
        
        if not hasattr(self, 'G') or not hasattr(self, 'Q'):
            self._build_influence_matrices()

        # Masks.
        free_surface = self.boundary.free_surface
        depth = self.boundary.depth
        bottom = self.boundary.depth
        cylinder = self.boundary.cylinder
        
        A = self.Q.copy().astype(np.complex128)
        q = np.zeros(self.boundary.number_of_elements).astype(np.complex128)
        
        A[:, free_surface] += -k * self.G[:, free_surface]
        A[:, depth] += 1j * k * self.G[:, depth]
        q[cylinder] = 1j * w * self.boundary.normals[cylinder, 1]  # Heave dof

        b = self.G @ q

        self.phi = np.linalg.solve(A, b)

    def get_radiation_coefficients_and_wave_amplitude(self, w):
        k = w**2/self.g  # wave number
        l = 2*np.pi/k  # wave length
        
        # Masks.
        free_surface = self.boundary.free_surface
        cylinder = self.boundary.cylinder

        normals = self.boundary.normals[cylinder]
        lengths = self.boundary.lengths[cylinder]

        nf = self.boundary.nf  # number of free surface elements on each side

        # Added mass, wave damping and wave amplitude.      
        self.solve(w)

        # Radiation coefficients.
        fy = 1j * self.rho * w * np.sum(self.phi[cylinder] * normals[:, 1] * lengths)
        az = -fy.real / w**2
        bz = fy.imag / w

        # Mask for free surface points at a distance larger than R + λ away from the cylinder.
        fsx = np.logical_and(np.abs(self.boundary.midpoints[:, 0]) > 1.0+l, free_surface)
        
        wsx = -1j * w / self.g * self.phi[fsx]

        # Wave amplitude.
        wa = np.abs(wsx).mean()

        return az, bz, wa