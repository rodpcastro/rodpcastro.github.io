import numpy as np

def get_domain_values(u, q, xv, yv, xb, yb, nx, ny, lm):
    """Find solution in any part of the domain."""

    n = len(xb) - 1
    s = np.zeros((len(xv), len(yv)))
    for i, x in enumerate(xv):
        for j, y in enumerate(yv):
            for k in range(n):
                [F, G] = findfg(x, y, xb[k], yb[k], nx[k], ny[k], lm[k])
                s[j, i] = s[j, i] + u[k]*G - q[k]*F

    return s
