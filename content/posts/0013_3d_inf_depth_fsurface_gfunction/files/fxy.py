import scipy as sc
import numpy as np
import mpmath as mp


mp.dps = 16
epsabs = 1e-10
epsrel = 1e-10


def fint(x, y):
    """F(X,Y) integral expression."""

    x2 = x**2
    ey = np.exp(-y)
    py = np.pi * ey
    dy = 2*ey

    phy0 = py*(sc.special.struve(0, x) + sc.special.y0(x))
    phy1 = py*(sc.special.struve(1, x) + sc.special.y1(x))
    phy2 = py*(sc.special.struve(2, x) + sc.special.yn(2, x))
    
    f = -dy * sc.integrate.quad(
        lambda t: np.exp(t) * (x2 + t**2)**-0.5,
        0.0,
        y,
        epsabs=epsabs,
        epsrel=epsrel,
    )[0] - phy0
    
    fx = dy*x * sc.integrate.quad(
        lambda t: np.exp(t) * (x2 + t**2)**-1.5,
        0.0,
        y,
        epsabs=epsabs,
        epsrel=epsrel,
    )[0] - dy + phy1
    
    fxx = dy * sc.integrate.quad(
        lambda t: np.exp(t) * (x2 + t**2)**-2.5 * (t**2 - 2*x2),
        0.0,
        y,
        epsabs=epsabs,
        epsrel=epsrel,
    )[0] + ey*x/3.0 + 0.5*(phy0 - phy2)
    
    return f, fx, fxx


def fsem(x, y):
    """F(X,Y) Series Expansion Method."""

    if x >= 9.5 and y > 0.5*x:
        # SEM4.
        if x < 14.0:
            nterms = 20
        else:
            nterms = 15
        
        return fsem4(x, y, nterms)

    elif x >= 6.5 and y <= 0.5*x:
        # SEM3.

        return fsem3(x, y, 13)

    elif y < 15.0 and y < 2*x:
        # SEM2.
        if y > 11.0:
            nterms = 42
        elif x > 7.0 or y > 8.0:
            nterms = 36
        elif x > 4.5 or y > 6.0:
            nterms = 30
        else:
            nterms = 26
        
        return fsem2(x, y, nterms)

    else:
        # SEM1.
        if 14.0 < y < 17.0:
            nterms = 19
        else:
            nterms = 15
        
        return fsem1(x, y, nterms)


def expei(x):
    """exp(-x) * Ei(x)."""

    if x > 40.0:
        sk = 1.0
        ex = 1.0
        for k in range(1, 24):
            sk *= k/x
            ex += sk
        ex *= 1/x
    else:
        ex = mp.exp(-x) * mp.ei(x)
    
    return ex


def fsem1(x, y, nterms=19):
    """Series Expansion Method 1."""

    eey = expei(y)

    xi = 1/x
    yi = 1/y
    
    qxi = 4*xi
    qx2 = 0.25*x*x
    
    tn1 = 1.0
    sn1 = 0.0
    sn2 = 0.0
    sn3 = 0.0
    for n in range(1, nterms+1):
        tn1 *= -qx2/(n*n)
        tn2 = tn1 * n
        tn3 = tn2 * (2*n-1)

        tm = yi
        sm = yi
        for m in range(2, 2*n+1):
            tm *= (m-1)*yi
            sm += tm
        sm -= eey
        
        sn1 += tn1 * sm
        sn2 += tn2 * sm
        sn3 += tn3 * sm

    f   = 2*(sn1 - eey)
    fx  = qxi * sn2
    fxx = qxi*xi * sn3

    return f, fx, fxx


def fsem2(x, y, nterms=42):
    """Series Expansion Method 2."""

    x2 = x**2
    y2 = y**2
    r2 = x2 + y2
    r = mp.sqrt(r2)
    xi = 1/x
    ri = 1/r

    rxi = r*xi
    rxi2 = rxi*xi
    r2xi2 = rxi2*r
    yxiri = y*xi*ri
    
    ey = mp.exp(-y)
    py = mp.pi * ey
    py0 = py * mp.bessely(0, x)
    py1 = py * mp.bessely(1, x)
    
    tn = 1.0
    sn1 = 1.0
    sn2 = 1.0
    sn3 = 0.0
    for n in range(1, nterms+1):
        tn *= x/n
        hg = mp.re(1j**-n * mp.hyp2f1(0.5, -0.5*n, 1.5, r2xi2))
        tg = tn * hg
        
        sn1 += (n+1) * tg
        sn2 += tg
        sn3 += n * tg

    f   = -py0 + 2*rxi2*(ey*sn1 - 1.0)
    fx  =  py1 + 2*(yxiri - rxi*ey*sn2)
    fxx =  py0 - py1*xi + 2*(yxiri*xi*(y2/r2 - 2.0 + y) - rxi2*ey*sn3)
    
    return f, fx, fxx


def fsem3(x, y, n3=13):
    """Series Expansion Method 3."""
    
    xi = 1/x
    xi2 = xi*xi
    xi3 = xi2*xi
    y2 = y*y
    yi = 1/y

    hxi2 = 0.5*xi2
    
    ey = mp.exp(-y)
    py = mp.pi * ey
    oy = 1.0 - ey
    
    phy0 = py * (mp.struveh(0, x) + mp.bessely(0, x))
    phy1 = py * (mp.struveh(1, x) + mp.bessely(1, x))

    tn = 1.0
    y2n = 1.0
    cn = oy
    sn1 = 0.0
    sn2 = 0.0
    sn3 = 0.0
    for n in range(1, n3+1):
        dn = 2*n
        tn *= -hxi2 * (dn-1) / n
        y2n *= y2
        cn = y2n*(1.0-dn*yi) + dn*(dn-1)*cn
        
        nc1 = tn * cn
        nc2 = (dn+1) * nc1
        nc3 = (dn+2) * nc2
        
        sn1 += nc1
        sn2 += nc2
        sn3 += nc3

    f   = -phy0 - 2*xi*(oy + sn1)
    fx  =  phy1 - 2*ey + 2*xi2*(oy + sn2)
    fxx =  phy0 - xi*phy1 - 2*xi3*(2*oy + sn3)
    
    return f, fx, fxx


def fsem4(x, y, n4=20):
    """Series Expansion Method 4."""

    x2 = x*x
    y2 = y*y
    r2 = x2 + y2
    r = mp.sqrt(r2)
    xi = 1/x
    yi = 1/y
    yi2 = yi*yi
    ri = 1/r
    ri2 = ri*ri
    ri3 = ri2*ri

    x2ri2 = x2*ri2
    y2ri2 = y2*ri2
    hyr = 0.5*y2ri2
    
    ey = mp.exp(-y)
    py = mp.pi * ey
    oy = 1 - ey
    eyi = ey*yi

    phy0 = py * (mp.struveh(0, x) + mp.bessely(0, x))
    phy1 = py * (mp.struveh(1, x) + mp.bessely(1, x))

    tn = -hyr
    b = 1.0
    b0 = oy*yi
    b1 = yi*((1.0 - 2*yi2)*ey - 2*(yi - yi2))
    sn1 = tn * b1
    sn2 = 3 * sn1
    sn3 = sn2 * (1 - 5*x2ri2)
    for n in range(2, n4+1):
        dn = 2*n
        tn *= -hyr * (dn-1) / n
        b *= -1.0
        bn = b*eyi + yi2*dn*((dn-1)*b1 + (dn-2)*b0)
        b0 = b1
        b1 = bn
        
        nb1 = tn * bn
        nb2 = nb1 * (dn+1)
        nb3 = nb2 * (1 - (dn+3)*x2ri2)
        
        sn1 += nb1
        sn2 += nb2
        sn3 += nb3

    f   = -phy0 - 2*ri*(oy + y*sn1)
    fx  =  phy1 -2*ey + 2*x*ri3*(oy + y*sn2)
    fxx =  phy0 - xi*phy1 + 2*ri3*(oy*(1 - 3*x2ri2) + y*sn3)
    
    return f, fx, fxx