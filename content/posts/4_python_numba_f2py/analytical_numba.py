# Extract from element.ipynb

import numpy as np
import numba

eps = np.finfo(np.float64).eps

class Element:
   
    @staticmethod
    @numba.jit('Tuple((f8, f8))(f8[:], f8)', nopython=True, cache=True)
    def _analytical_numba(field_local, element_length):
        """Get influence coefficients at a field point using analytical integration."""

        a = 0.5 * element_length
        x, y = field_local

        if np.abs(y) <= element_length * eps \
           and np.abs(np.abs(x) - a) <= element_length * eps:
            G = a / np.pi * (np.log(2*a) - 1)
            Q = 0.0
        else:
            xpa = x + a
            xma = x - a

            r1 = np.sqrt(xma**2 + y**2)
            r2 = np.sqrt(xpa**2 + y**2)
            t1 = np.arctan2(y, xma)
            t2 = np.arctan2(y, xpa)

            G = 0.5 / np.pi * (
                y * (t1 - t2) - xma * np.log(r1) + xpa * np.log(r2) - 2 * a
            )

            if np.abs(y) <= element_length * eps:
                # Q is discontinuous in |x| < a and y = 0.
                Q = 0.0
            else:
                Q = -0.5 / np.pi * (t1 - t2)

        return G, Q
