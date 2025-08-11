import numpy as np
from scipy.integrate import quad

def intf(t, xi, eta, xk, yk, nkx, nky, lk):
    """Integrand of F."""

    return np.log((xk - t*lk*nky - xi)**2 + (yk + t*lk*nkx - eta)**2)

def intg(t, xi, eta, xk, yk, nkx, nky, lk):
    """Integrand of G."""

    return (nkx*(xk - t*lk*nky - xi) + nky*(yk + t*lk*nkx - eta)) / ((xk - t*lk*nky - xi)**2 + (yk + t*lk*nkx - eta)**2)

def findfg(xi, eta, xk, yk, nkx, nky, lk):
    """Functions F and G"""

    F = (lk/(4.0*np.pi)) * quad(lambda t: intf(t, xi, eta, xk, yk, nkx, nky, lk), 0, 1, epsabs=1e-8)[0]
    G = (lk/(2.0*np.pi)) * quad(lambda t: intg(t, xi, eta, xk, yk, nkx, nky, lk), 0, 1, epsabs=1e-8)[0]

    return F, G
