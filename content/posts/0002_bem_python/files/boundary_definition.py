import numpy as np

def define_boundary(nl):
    """Define 1x1 Square region uniformly discretized in 4*nl elements.

    Parameters
    ----------
    nl : int
        Number of elements per side.

    Returns
    -------
    xb : numpy.ndarray[float]
        X coordinate of boundary node points.
    yb : numpy.ndarray[float]
        Y coordinate of boundary node points.
    bt : numpy.ndarray[int]
        Boundary condition type by element.
    bv : numpy.ndarray[float]
        Boundary condition value by element.
    """

    # y = 0 and 0 ≤ x ≤ 1.
    xb1 = np.linspace(0, 1, nl+1)[:-1]
    yb1 = np.zeros(nl+1)[:-1]
    bt1 = np.ones(nl, dtype=np.int8)
    bv1 = np.zeros(nl)
    
    # x = 1 and 0 ≤ y ≤ 1.
    xb2 = np.ones(nl+1)[:-1]
    yb2 = np.linspace(0, 1, nl+1)[:-1]
    bt2 = np.zeros(nl, dtype=np.int8)
    ym = np.linspace(0.5/nl, 1-0.5/nl, nl)
    bv2 = np.cos(np.pi*ym)

    # y = 1 and 0 ≤ x ≤ 1.
    xb3 = np.linspace(1, 0, nl+1)[:-1]
    yb3 = np.ones(nl+1)[:-1]
    bt3 = np.ones(nl, dtype=np.int8)
    bv3 = np.zeros(nl)
    
    # x = 0 and 0 ≤ y ≤ 1.
    xb4 = np.zeros(nl+1)
    yb4 = np.linspace(1, 0, nl+1)
    bt4 = np.zeros(nl, dtype=np.int8)
    bv4 = np.zeros(nl)

    xb = np.concatenate((xb1, xb2, xb3, xb4))
    yb = np.concatenate((yb1, yb2, yb3, yb4))
    bt = np.concatenate((bt1, bt2, bt3, bt4))
    bv = np.concatenate((bv1, bv2, bv3, bv4))
    
    return xb, yb, bt, bv
