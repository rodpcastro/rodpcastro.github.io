---
date: '2025-07-19'
draft: true
title: 'Using Julia to compute the analytical solution of a two-dimensional potential flow with free-surface'
author: 'Rodrigo Castro'
summary: 'Hydrodynamic coefficients and the wave produced by a circular cylinder oscilating in heave are obtained by evaluation of analytical expressions. Experimental data is also presented for validation.'
tags: ['Julia', 'Potential flow']
---

## Introduction
This post showcases my first steps in learning the Julia programming language. The chosen application is the evaluation of the analytical expressions describing the two-dimensional potential flow of a circular cylinder oscilating vertically on the free surface of a fluid.

The Julia code is not presented here, but left as an [appendix](#appendices).

## Methods
The potential flow studied here consists of a cylinder of circular section with its axis parallel to the free surface. If the cylinder is given a forced simple harmonic motion of small amplitude about its initial position, waves travel away from it. At a distance of a few wave-lengths from the cylinder, the motion is described by a singular regular wave-train of amplitude proportional to the amplitude of the forced oscillation. From the potential or the stream function, it's possible to deduce the wave amplitude at a distance from the cylinder and the added mass and damping due to the fluid motion.

### Analytical expressions
The analytical expression here presented were derived by Ursell (1949). The following figure depicts a cylinder of radius $a$ and the coordinate systems used. 

{{< figure src="images/cylinder.svg" alt="Cylinder" align="center" >}}

Since the motion is symmetric about the y-axis, it is sufficient to consider the quadrant $0 \le \theta \le \frac{\pi}{2}$. The potential $\phi$ and stream $\psi$ functions satisfy Laplace's equation:

$$\eq{
\frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2} = 0,
}$$

$$\eq{
\frac{\partial^2 \psi}{\partial x^2} + \frac{\partial^2 \psi}{\partial y^2} = 0.
}$$

The boundary condition of the free surface is describe by

$$\eq{
K\phi + \frac{\partial \phi}{\partial y} = 0, \quad \theta=\frac{\pi}{2}, \quad r > a,
}$$

where $K = \sigma^2 / g$ is the wavenumber. Also, by symmetry,

$$\eq{
\frac{\partial \phi}{\partial \theta} = 0, \quad \theta=0.
}$$

The boundary condition on the cylinder is expressed by

$$\eq{
\psi = l \sigma a \sin(\sigma t + \epsilon)\sin\theta \quad \text{on} \quad r = a,
}$$

where $l$ is the amplitude of the forced vertical motion and $\sigma$ is the frequency.

The required solution must satisfy the boundary conditions $(3)$, $(4)$ and $(5)$ and the radiation condition at infinity. The solution derived by Ursell (1949) is composed of a source potential in the origin and a sum of multipole potentials, for which the coefficients are chosen in such a way that all boundary conditions are satisfied.

The potential function $\phi$, solution of the described boundary value problem is expressed in the form

$$\eq{\begin{split}
\frac{\pi \sigma \phi}{gb} = 
\cos(\sigma t) \left( \Phi_c(K, r, \theta) + \sum_{m=1}^{\infty} p_{2m} \Phi_m(K, r, \theta) \right) + \\
\sin(\sigma t) \left( \Phi_s(K, r, \theta) + \sum_{m=1}^{\infty} q_{2m} \Phi_m(K, r, \theta) \right), \phantom{,} 
\end{split}}$$

where

$$\eq{
\Phi_c(K, r, \theta) = \pi e^{-K r \cos\theta} \cos(K r \sin\theta),
}$$

$$\eq{\begin{split}
\Phi_s(K, r, \theta) = -\int_{0}^{\infty} \frac{e^{-k r \sin\theta}}{K^2+k^2} \left[ k\cos(k r \cos\theta) - K\sin(k r \cos\theta) \right]dk + \\
\quad + \pi e^{-K r \cos\theta} \sin(K r \sin\theta),
\end{split}}$$

$$\eq{
\Phi_m(K, r, \theta) = \left(\frac{a}{r}\right)^{2m} \left[\cos(2m\theta) + \frac{K r}{2m-1}\cos(2m\theta - \theta) \right].
}$$

Similarly, the stream function $\psi$ is

$$\eq{\begin{split}
\frac{\pi \sigma \psi}{gb} = 
\cos(\sigma t) \left( \Psi_c(K, r, \theta) + \sum_{m=1}^{\infty} p_{2m} \Psi_m(K, r, \theta) \right) + \\
\sin(\sigma t) \left( \Psi_s(K, r, \theta) + \sum_{m=1}^{\infty} q_{2m} \Psi_m(K, r, \theta) \right), \phantom{,} 
\end{split}}$$

where

$$\eq{
\Psi_c(K, r, \theta) = \pi e^{-K r \cos\theta} \sin(K r \sin\theta),
}$$

$$\eq{\begin{split}
\Psi_s(K, r, \theta) = \int_{0}^{\infty} \frac{e^{-k r \sin\theta}}{K^2+k^2} \left[ k\sin(k r \cos\theta) + K\cos(k r \cos\theta) \right]dk - \\
\quad - \pi e^{-K r \cos\theta} \cos(K r \sin\theta),
\end{split}}$$

$$\eq{
\Psi_m(K, r, \theta) = \left(\frac{a}{r}\right)^{2m} \left[\sin(2m\theta) + \frac{K r}{2m-1}\sin(2m\theta - \theta) \right].
}$$

The coefficients $p_{2m}$ and $q_{2m}$ are found by solving $(6)$ or $(10)$ for the boundary condition $(5)$. The expressions for finding these coefficients are given below using expression $(10)$ as starting point.

$$\eq{
\sum_{m=1}^{\infty} p_{2m} f_{2m}(\theta) = \Psi_c(K, a, \theta) - \Psi_c(K, a, \frac{\pi}{2}) \sin\theta,
}$$

$$\eq{
\sum_{m=1}^{\infty} q_{2m} f_{2m}(\theta) = \Psi_s(K, a, \theta) - \Psi_s(K, a, \frac{\pi}{2}) \sin\theta.
}$$

Following the procedure used by Ursell, each of the two equations above is written as an overdetermined system, solved by the least squares method. The system is built by choosing 10 values of $\theta$ ranging from $0$ to $\frac{pi}{2}$, and 6 values of $m$ ranging from $1$ to $6$.

After that, one can define

$$\eq{
A = \Psi_c(K, a, \frac{\pi}{2}) + \sum_{m=1}^{\infty} \frac{(-1)^{m-1} K a}{2m-1} p_{2m}, 
}$$

$$\eq{
B = \Psi_s(K, a, \frac{\pi}{2}) + \sum_{m=1}^{\infty} \frac{(-1)^{m-1} K a}{2m-1} q_{2m}.
}$$

Then the stream function on the cylinder is given by

$$\eq{
\psi = \frac{g b}{\pi \sigma}\left[ A \cos(\sigma t) + B \sin(\sigma t) \right].
}$$

Combining the equation above with the boundary condition defined in $(5)$, the ratio between the radiated wave amplitude $b$ and the forced oscilation amplitude $l$ is

$$\eq{
\frac{b}{l} = \frac{\pi K a}{\sqrt{A^2 + B^2}}.
}$$

Omitting further calculation details, the dimensionless added mass and wave damping coefficients are defined by

$$\eq{
m_z = \frac{4}{\pi} \frac{M_0 B + N_0 A}{A^2 + B^2},
}$$

$$\eq{
b_z = \frac{4 \sigma \sqrt{a}}{\pi \sqrt{g}} \frac{M_0 A - N_0 B}{A^2 + B^2},
}$$

where

$$\eq{
M_0 = \int_{0}^{\frac{\pi}{2}} \Phi_s(K, a, \theta) \cos\theta\,d\theta + \frac{\pi}{4} K a q_2 + \sum_{m=1}^{\infty} \frac{(-1)^{m-1}q_{2m}}{4m^2 - 1},
}$$

$$\eq{
N_0 = \int_{0}^{\frac{\pi}{2}} \Phi_c(K, a, \theta) \cos\theta\,d\theta + \frac{\pi}{4} K a p_2 + \sum_{m=1}^{\infty} \frac{(-1)^{m-1}p_{2m}}{4m^2 - 1}.
}$$

## Results
The results obtained from the evaluation of the analytical expressions are compared with experimental data obtained by Vugts (1968). In the following plots, $\omega = \sigma\sqrt{\frac{a}{g}}$ is the dimensionless frequency.

In general, analytical expressions and experimental data match well, except the added mass at low frequencies which, by the analytical formula, tend to infinity as the frequency tends to zero.

{{< figure src="images/rao.svg" alt="Wave/Heave amplitudes ratio" align="center" >}}

{{< figure src="images/azz.svg" alt="Heave added mass" align="center" >}}

{{< figure src="images/bzz.svg" alt="Heave wave damping" align="center" >}}

## Conclusion
This small project was a good start for learning basic elements of the Julia programming language. Most importantly, the results here presented will be reused in future posts. Next time, analytical and experimental data will be compared with results obtained by the boundary element method.

## References
1. F. Ursell. 1949. On the heaving motion of a circular cylinder on the surface of a fluid. The Quarterly Journal of Mechanics and Applied Mathematics, 2, 2 (1949), 218–231. https://doi.org/10.1093/qjmam/2.2.218
2. J. H. Vugts. 1968. The hydrodynamic coefficients for swaying, heaving and rolling cylinders in a free surface. International Shipbuilding Progress, 15, 167 (1968), 251–276. https://doi.org/10.3233/ISP-1968-1516702

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[Julia]: https://julialang.org/
