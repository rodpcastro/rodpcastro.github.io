import numpy as np

def solve_bem(xb, yb, bt, bv, xm, ym, lm, nx, ny):
    """Build and solve system of equations.

    Returns
    -------
    u : numpy.ndarray[float]
        Function at the boundaries.
    q : numpy.ndarray[float]
        Derivative at the boundaries.
    """

    n = len(xb) - 1
    b = np.zeros(n)
    A = np.zeros((n,n))
    for m in range(n):
        for k in range(n):
            if k == m:
                G = 0.0
                F = lm[k]/(2.0*np.pi) * (np.log(lm[k]/2.0) - 1.0)
                d = 1.0
            else:
                F, G = findfg(xm[m], ym[m], xb[k], yb[k], nx[k], ny[k], lm[k])
                d = 0.0
            if bt[k] == 0:
                A[m, k] = -F
                b[m] = b[m] + bv[k]*(-G + 0.5*d)
            else:
                A[m, k] = G - 0.5*d
                b[m] = b[m] + bv[k]*F
    
    z = np.linalg.solve(A, b)
    
    # Assign approximate boundary values accordingly.
    u, q = [np.empty(n) for _ in range(2)]
    for m in range(n):
        u[m] = (1 - bt[m])*bv[m] + bt[m]*z[m]
        q[m] = (1 - bt[m])*z[m] + bt[m]*bv[m]

    return u, q
