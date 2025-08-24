import numpy as np
from scipy.special import exp1
from twodubem.green import Green

class FreeSurface(Green):
    """Infinite-depth free-surface Green Function.

    Parameters
    ----------
    w : float
        Wave frequency.
    body : Polygon, default=None
        Body to compute influence matrices.
    """
    
    def __init__(self, w, body=None):
        self.w = w
        self.g = 9.81
        self.K = w**2 / self.g
        if body is not None:
            self._build_influence_matrices(body)
    
    @staticmethod
    def eval(field_point, source_point, K):
        """Infinite-depth free-surface Green function."""

        sx, sz = source_point
        fx, fz = field_point

        x1 = fx - sx
        z1 = fz - sz
        z3 = fz + sz

        R = np.abs(x1)
        v1 = np.abs(z1)
        v3 = np.abs(z3)

        # Auxilary variables d.
        d1 = R**2
        d2 = v1**2
        d3 = v3**2
        d4 = d1 + d2
        d5  = d4**2
        d6 = d1 + d3
        d7 = d6**2
        d8 = 2*R

        r1 = np.sqrt(d4)
        r3 = np.sqrt(d6)

        X = K*R
        V3 = K*v3
        Z = V3 - 1j*X

        # Auxilary variables k.
        k1 = 2*K
        k2 = k1*K

        # Auxilary variables e.
        e1 = np.exp(-Z)
        e2 = e1*exp1(-Z)
        e3 = e2 + 1/Z
        e4 = e3 + 1/Z**2
        e5 = 2*np.pi*e1
        e6 = K*e5
        e7 = K*e6

        # Auxilary variables s.
        sx1 = np.sign(x1)
        sz1 = np.sign(z1)

        if X <= 1:
            # Near field
            G = np.log(K*r1) + np.log(K*r3) - 2 * (e2.real + np.log(np.abs(Z))) - 1j*e5
        else:
            # Far field
            G = np.log(r1/r3) - 2*e2.real - 1j*e5

        Gx = sx1 * (R/d4 - R/d6 + k1*e3.imag + e6)
        Gz = sz1 * v1/d4 + v3/d6 - k1*e3.real - 1j*e6
        gradG = np.array([Gx, Gz], dtype=np.complex128)

        Gxx = (d2 - d1)/d5 + (d1 - d3)/d7 + k2*e4.real + 1j*e7
        Gzz = -Gxx
        Gxz = -sx1 * (sz1 * v1*d8/d5 + v3*d8/d7 - k2*e4.imag - e7)
        hessG = np.array([[Gxx, Gxz], [Gxz, Gzz]], dtype=np.complex128)

        return G, gradG, hessG

    def get_element_influence_coefficients(self, element, point):
        a = 0.5 * element.length

        # 4-point Gauss-Legendre quadrature roots and weights.
        roots = [0.3399810435848563, 0.8611363115940526]
        weights = [0.6521451548625461, 0.3478548451374538]

        G = 0.0 + 1j*0.0
        gradG = np.zeros(2, dtype=np.complex128)
        hessG = np.zeros((2, 2), dtype=np.complex128)
        for i in range(len(roots)):
            element_point_p = element.get_point_global_coordinates(
                np.array([a * roots[i], 0.0])
            )
            element_point_m = element.get_point_global_coordinates(
                np.array([-a * roots[i], 0.0])
            )

            g_p, gradg_p, hessg_p = self.eval(element_point_p, point, self.K)
            g_m, gradg_m, hessg_m = self.eval(element_point_m, point, self.K)

            G += weights[i] * (g_p + g_m)
            gradG += weights[i] * (gradg_p + gradg_m)
            hessG += weights[i] * (hessg_p + hessg_m)

        G = a * G
        Q = a * gradG @ element.normal
        gradG = -a * gradG
        gradQ = -a * hessG @ element.normal

        return G, Q, gradG, gradQ

    def _build_influence_matrices(self, body):
        n = body.number_of_elements
        self.G = np.empty((n, n), dtype=np.complex128)
        self.Q = np.empty((n, n), dtype=np.complex128)

        for i, source_element in enumerate(body.elements):
            source_point = source_element.node
            for j, field_element in enumerate(body.elements):
                g, q, _, _ = self.get_element_influence_coefficients(
                        field_element,
                        source_point,
                )
                self.G[i, j] = g
                self.Q[i, j] = q
                if i == j:
                    if i < body.number_of_body_elements:
                        self.Q[i, j] += -np.pi
                    else:
                        self.Q[i, j] += 2*np.pi