import numpy as np
from twodubem.solver import Solver


class WaveSolver(Solver):

    def __init__(self, body, green):
        self.body = body
        self.green = green
        self.g = green.g
        self.w = green.w
        self.K = green.K
        self.L = 2*np.pi / self.K
        self.rho = 1.0
        self._build_boundary_condition_vector()

    def _build_boundary_condition_vector(self):
        self.q = np.zeros(self.body.number_of_elements, dtype=np.complex128)

    def solve(self):
        b = self.green.G @ self.q
        self.phi = np.linalg.solve(self.green.Q, b)

    def get_solution(self, X, Z):
        """Get solution for array of points."""

        n = self.body.number_of_elements
        x = X.ravel()
        z = Z.ravel()
        u = np.empty(x.shape, dtype=np.complex128)
        v = np.empty((*x.shape, 2), dtype=np.complex128)
        G = np.empty(n, dtype=np.complex128)
        Q = np.empty(n, dtype=np.complex128)
        gradG = np.empty((2, n), dtype=np.complex128)
        gradQ = np.empty((2, n), dtype=np.complex128)
        dpi = 2*np.pi

        for i in range(len(x)):

            field_point = np.array([x[i], z[i]])

            for j, source_element in enumerate(self.body.elements):
                g, q, gradg, gradq = self.green.get_element_influence_coefficients(
                    source_element,
                    field_point,
                )
                G[j] = g
                Q[j] = q
                gradG[:, j] = gradg
                gradQ[:, j] = gradq
            
            u[i] = Q @ self.phi - G @ self.q
            v[i] = gradQ @ self.phi - gradG @ self.q

        phi = u.reshape(X.shape) / dpi
        gradphi = v.reshape((*X.shape, 2)) / dpi

        return phi, gradphi

    def __add__(self, other):
        if not isinstance(other, WaveSolver):
            raise ValueError("Operands must be WaveSolver")
        elif other.body is not self.body:
            raise ValueError("Operands must have the same boundary")
        elif other.green is not self.green:
            raise ValueError("Operands must have the same Green function")
            
        solver_sum = WaveSolver(self.body, self.green)
        solver_sum.phi = self.phi + other.phi
        solver_sum.q = self.q + other.q

        return solver_sum


class RadiationSolver(WaveSolver):
    """Solver for the radiation problem."""

    def _build_boundary_condition_vector(self):
        nd = len(self.body.dofs)
        ne = self.body.number_of_elements
        body = self.body._mask_body
        
        self.q = np.zeros((ne, nd), dtype=np.complex128)
        
        for i, dof in enumerate(self.body.dofs.keys()):
            self.q[body, i] = -1j * self.w * self.body.dofs[dof]

    def compute_radiation_coefficients(self):
        nd = len(self.body.dofs)
        body = self.body._mask_body
        lengths = self.body.lengths[body]
        f = np.zeros((nd, nd), dtype=np.complex128)

        for i, dofi in enumerate(self.body.dofs):
            for j, dofj in enumerate(self.body.dofs):
                f[i, j] = 1j * self.w * self.rho * np.sum(self.phi[body, j] * self.body.dofs[dofi] * lengths)

        self.added_mass = f.real / self.w**2
        self.radiation_damping = f.imag / self.w

    def get_solution(self, X, Z):
        nd = len(self.body.dofs)
        ne = self.body.number_of_elements
        x = X.ravel()
        z = Z.ravel()
        u = np.empty((*x.shape, nd), dtype=np.complex128)
        v = np.empty((*x.shape, 2, nd), dtype=np.complex128)
        G = np.empty(ne, dtype=np.complex128)
        Q = np.empty(ne, dtype=np.complex128)
        gradG = np.empty((2, ne), dtype=np.complex128)
        gradQ = np.empty((2, ne), dtype=np.complex128)
        dpi = 2*np.pi

        for i in range(len(x)):

            field_point = np.array([x[i], z[i]])

            for j, source_element in enumerate(self.body.elements):
                g, q, gradg, gradq = self.green.get_element_influence_coefficients(
                    source_element,
                    field_point,
                )
                G[j] = g
                Q[j] = q
                gradG[:, j] = gradg
                gradQ[:, j] = gradq
            
            u[i] = Q @ self.phi - G @ self.q
            v[i] = gradQ @ self.phi - gradG @ self.q

        phi = u.reshape((*X.shape, nd)) / dpi
        gradphi = v.reshape((*X.shape, 2, nd)) / dpi

        return phi, gradphi


class DiffractionSolver(WaveSolver):
    """Solver for the diffraction problem."""

    def _build_boundary_condition_vector(self):
        ne = self.body.number_of_elements
        body = self.body._mask_body
        normals = self.body.normals[body]

        X = self.body.midpoints[body, 0]
        Z = self.body.midpoints[body, 1]
        
        phi0, gradphi0 = self.incident_wave_potential(X, Z, self.w, self.g)

        self.q = np.zeros(ne, dtype=np.complex128)
        self.q[body] = -np.sum(gradphi0 * normals, axis=1)
        
        b = self.green.G @ self.q

        self.phi0 = phi0
        self.phi = np.linalg.solve(self.green.Q, b)

    def compute_exciting_forces(self):
        nd = len(self.body.dofs)
        body = self.body._mask_body
        normals = self.body.normals[body]
        lengths = self.body.lengths[body]
        phi = self.phi0 + self.phi[body]
        f = np.zeros(nd, dtype=np.complex128)

        for i, dof in enumerate(self.body.dofs):
            f[i] = 1j * self.w * self.rho * np.sum(phi * self.body.dofs[dof] * lengths)

        self.force = f

    @staticmethod
    def incident_wave_potential(X, Z, w, g=9.81):
        K = w**2/g
        exz = np.exp(K*(Z + 1j*X))
        phi = 1j*g/w * exz
        phix = -w * exz
        phiz = -1j * phix
        gradphi = np.stack((phix, phiz), axis=-1)

        return phi, gradphi