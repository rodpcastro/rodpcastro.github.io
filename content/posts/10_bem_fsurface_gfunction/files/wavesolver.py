import numpy as np
from twodubem.solver import Solver
from wavegreen import FreeSurfaceGreenFunction


class WaveSolver(Solver):
    """Solver for the radiation and diffraction problems."""
    
    def __init__(self, body, w):
        self.boundary = body
        self.green = FreeSurfaceGreenFunction()
        self.g = 9.81  # Acceleration of gravity
        self.rho = 1.0  # Water density
        self.w = w
        self.K = self.w**2 / self.g
        self.L = 2*np.pi / self.K

    def _build_influence_matrices(self):
        """Build influence coefficients matrices."""
 
        n = self.boundary.number_of_elements
        self.G = np.empty((n, n), dtype=np.complex128)
        self.Q = np.empty((n, n), dtype=np.complex128)

        for i, source_element in enumerate(self.boundary.elements):
            source_point = source_element.node
            for j, field_element in enumerate(self.boundary.elements):
                g, q, _, _ = self.green.get_line_element_influence_coefficients(
                        field_element,
                        source_point,
                        self.K,
                )
                self.G[i, j] = g
                self.Q[i, j] = q
                if i == j:
                    self.Q[i, j] += -np.pi

    def solve_radiation_problem(self):
        """Solve the radiation potential of oscilation in heave."""

        if not hasattr(self, 'G') or not hasattr(self, 'Q'):
            self._build_influence_matrices()

        # Heave dof.
        self.qr = -1j * self.w * self.boundary.normals[:, 1]

        b = self.G @ self.qr

        self.phi_radiation = np.linalg.solve(self.Q, b)

    def solve_diffraction_problem(self):
        """Solve the diffraction potential."""

        if not hasattr(self, 'G') or not hasattr(self, 'Q'):
            self._build_influence_matrices()

        phi, gradphi = self.incident_wave_potential(self.boundary.midpoints[:, 0],
                                                    self.boundary.midpoints[:, 1],
                                                    self.w,
                                                    self.g)

        self.qd = -np.sum(gradphi * self.boundary.normals, axis=1)
        
        b = self.G @ self.qd

        self.phi_incident = phi
        self.phi_diffraction = np.linalg.solve(self.Q, b)
    
    def get_radiation_coefficients(self):
        normals = self.boundary.normals
        lengths = self.boundary.lengths

        fz = 1j * self.w * self.rho * np.sum(self.phi_radiation * normals[:, 1] * lengths)
        az = fz.real / self.w**2
        bz = fz.imag / self.w

        return az, bz

    def get_potentials(self, X, Z):
        n = self.boundary.number_of_elements
        x = X.ravel()
        z = Z.ravel()
        zr = np.empty(x.shape, dtype=np.complex128)
        zd = np.empty(x.shape, dtype=np.complex128)
        wr = np.empty((*x.shape, 2), dtype=np.complex128)
        wd = np.empty((*x.shape, 2), dtype=np.complex128)
        G = np.empty(n, dtype=np.complex128)
        Q = np.empty(n, dtype=np.complex128)
        gradG = np.empty((2, n), dtype=np.complex128)
        gradQ = np.empty((2, n), dtype=np.complex128)
        dpi = 2*np.pi

        for i in range(len(x)):

            field_point = np.array([x[i], z[i]])
            
            for j, source_element in enumerate(self.boundary.elements):
                g, q, gradg, gradq = self.green.get_line_element_influence_coefficients(
                        source_element,
                        field_point,
                        self.K,
                )
                G[j] = g
                Q[j] = q
                gradG[:, j] = gradg
                gradQ[:, j] = gradq
            
            zr[i] = Q @ self.phi_radiation - G @ self.qr
            zd[i] = Q @ self.phi_diffraction - G @ self.qd
    
            wr[i] = gradQ @ self.phi_radiation - gradG @ self.qr
            wd[i] = gradQ @ self.phi_diffraction - gradG @ self.qd

        phir = zr.reshape(X.shape) / dpi
        phid = zd.reshape(X.shape) / dpi
        
        gradphir = wr.reshape((*X.shape, 2)) / dpi
        gradphid = wd.reshape((*X.shape, 2)) / dpi
        
        return phir, phid, gradphir, gradphid
    
    def get_wave_amplitude(self):
        phir, _, _, _ = self.get_potentials(np.array([2.0*self.L]), np.array([0.0]))

        zeta = np.abs(1j * phir[0]) * self.w / self.g

        return zeta

    def get_forces(self):
        normals = self.boundary.normals
        lengths = self.boundary.lengths
        phi_total = self.phi_incident + self.phi_diffraction

        fz = 1j * self.w * self.rho * np.sum(phi_total * normals[:, 1] * lengths)

        return fz

    @staticmethod
    def incident_wave_potential(X, Z, w, g=9.81):
        k = w**2/g
        exz = np.exp(k*(Z + 1j*X))
        phi = 1j*g/w * exz
        phix = -w * exz
        phiz = -1j * phix
        gradphi = np.stack((phix, phiz), axis=-1)

        return phi, gradphi