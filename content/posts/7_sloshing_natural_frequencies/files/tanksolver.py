import numpy as np
from twodubem.solver import Solver
from twodubem.laplace import Laplace


class SloshingSolver(Solver):
    """Solver for the linear sloshing problem.

    Parameters
    ----------
    tank : Polygon
        Boundary that represents the tank.
    """

    def __init__(self, tank):
        self.boundary = tank
        self.green = Laplace()
        self.method = 'constant'
        self.g = 9.81

    def solve_eigenvalue_problem(self):
        """Solve the eigenvalue problem.
        
        The solution of the eigenvalue problem returns natural sloshing frequencies and modes.
        """

        self._build_influence_matrices()
        
        n = self.boundary.number_of_elements
        nf = self.boundary.number_of_free_surface_elements
        nr = n - nf
        A = np.empty((n, n), dtype=np.float64)
        B = np.empty((n, n), dtype=np.float64)

        A[:, :nr] =  self.Q[:, :nr]
        A[:, nr:] = -self.G[:, nr:]
        B[:, :nr] =  self.G[:, :nr]
        B[:, nr:] = -self.Q[:, nr:]

        C = np.zeros_like(B)
        C[nr:, nr:] = np.eye(nf)

        R = B @ C

        T = np.linalg.inv(A) @ R
        T1 = T[nr:, nr:]

        s, q = np.linalg.eig(T1)

        sorting_array = np.argsort(s)
        s_sorted = s[sorting_array]
        q_sorted = q[:, sorting_array]
        
        self.natural_frequencies = np.sqrt(self.g * s_sorted[1:])
        self.natural_modes = q_sorted[:, 1:]