import numpy as np

def get_domain_values(u, q, xv, yv, xb, yb, nx, ny, lm):
    """Find solution in any part of the domain."""

    n = len(xb) - 1
    s = np.zeros((len(yv), len(xv)))  # To match numpy.meshgrid's default indexing.
    for j, y in enumerate(yv):
        for i, x in enumerate(xv):
            for k in range(n):
                [F, G] = findfg(x, y, xb[k], yb[k], nx[k], ny[k], lm[k])
                s[j, i] = s[j, i] + u[k]*G - q[k]*F

    return s
