import numpy as np

def elements_properties(xb, yb):
    """Find elements midpoints, lengths and unit normal vectors.

    Returns
    -------
    xm : numpy.ndarray[float]
        X coordinate of elements midpoints.
    ym : numpy.ndarray[float]
        Y coordinate of elements midpoints.
    lm : numpy.ndarray[float]
        Length of elements.
    nx : numpy.ndarray[float]
        X component of element normal vector.
    ny : numpy.ndarray[float]
        Y component of element normal vector.
    """

    n = len(xb) - 1
    xm, ym, lm, nx, ny = [np.empty(n) for _ in range(5)]
    for i in range(n):
        xm[i] = 0.5*(xb[i] + xb[i + 1])
        ym[i] = 0.5*(yb[i] + yb[i + 1])
        lm[i] = np.sqrt((xb[i+1] - xb[i])**2 + (yb[i + 1] - yb[i])**2)
        nx[i] = (yb[i+1] - yb[i]) / lm[i]
        ny[i] = (xb[i] - xb[i+1]) / lm[i]

    return xm, ym, lm, nx, ny