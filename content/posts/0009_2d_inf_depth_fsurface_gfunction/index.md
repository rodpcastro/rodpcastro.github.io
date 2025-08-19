---
date: '2025-08-08'
draft: false
title: 'The two-dimensional infinite-depth free-surface Green function'
author: 'Rodrigo Castro'
summary: 'Expressions for the infinite-depth Green function for diffraction and radiation of regular waves by two-dimensional structures.'
tags: ['Green Function', 'Potential Flow']
---

## Introduction
This post is dedicated to presenting the mathematical expressions that define the two-dimensional infinite-depth free-surface Green function. This Green function will be used in future posts to solve the diffraction and radiation problems of a two-dimensional structure in regular waves.

## Expressions
First, the quantities are defined and the Green function is presented. Then, expressions aimed for computational implementation are introduced.

### Free-surface Green function
The two-dimensional coordinate system, shown in the figure below, is defined with the $x$-axis horizontal and coincident with the mean free surface, the $z$-axis perpendicular to $x$, with $z = 0$ on the mean free surface and depth increasing with the $-z$ direction.

{{< figure src="images/coordinate_system.svg" alt="Coordinate system" align="center" >}}

The Green function $G(\mathrm{P}, \mathrm{Q})$ describes the spacial component of a velocity potential of the form $\Phi(x, z, t) = \mathfrak{Re}[G(\mathrm{P}, \mathrm{Q}) e^{-\mathrm{i} \omega t}]$ induced at a field point $\mathrm{Q} = (x, z)$  by a pulsating source located at $\mathrm{P} = (\xi, \zeta)$, oscillating at angular frequency $\omega$, in waters of infinite depth, . The Green function satisfies the Laplace equation in the fluid domain, the linearized free-surface boundary condition, and a radiation condition in the far field, given by

$$\eq{
\nabla^2 G = \frac{1}{2\pi} \delta(x-\xi) \delta(z-\zeta), \quad \text{in the fluid},
}$$

$$\eq{
\frac{\partial G}{\partial z} = K G, \quad \text{on } z = 0,
}$$

$$\eq{
\frac{\partial G}{\partial x} = \pm \mathrm{i} K G, \quad \text{for } K (x-\xi) \to \pm \infty,
}$$

where $\delta(.)$ is the Dirac delta function, $K = \omega^2 / g$ is the infinite-depth wave number, and $g$ is the acceleration of gravity. The following quantities are defined

$$\eq{
v_1 = |z - \zeta|, \quad v_3 = |z + \zeta|,
}$$

$$\eq{
R = |x - \xi|, \quad r_1 = \sqrt{R^2 + v_1^2}, \quad r_3 = \sqrt{R^2 + v_3^2},
}$$

which are depicted in the previous image. $r_1$ is the distance from $\mathrm{Q}$ to $\mathrm{P}$, $r_3$ is the distance from $\mathrm{Q}$ to $\mathrm{P}''$, the image of $\mathrm{P}$ relative to the mean free surface. We also define the dimensionless quantities

$$\eq{
X = K R, \quad V_3 = K v_3, \quad Z = V_3 - \mathrm{i} X.
}$$

Then, the free-surface Green function for infinite-depth is given by

$$\eq{\begin{split}
G^{\text{deep}} &= \log(K r_1) + \log(K r_3) \\
&\quad - 2 \left\{ \mathfrak{Re}\left[ e^{-Z} \mathrm{E}_1(-Z) \right] + \log(|Z|) \right\} - 2 \pi \mathrm{i} e^{-Z}
\end{split}}$$

For $X \gt 1$, the following simplification is suggested

$$\eq{
G^{\text{deep}} = \log\left( \frac{r_1}{r_3} \right) - 2 \mathfrak{Re}\left[ e^{-Z} \mathrm{E}_1(-Z) \right] - 2 \pi \mathrm{i} e^{-Z}
}$$

### Expressions for computational use
Now we introduce the Green function and its derivatives of first and second order calculated as a function of several auxilary variables, in a way that is useful for computational evaluation. All derivatives were evaluated manually and validated with [SymPy]. In the following definitions, auxilary quantities are given in unnumbered equations.

$$\begin{split}
x_1 &= x - \xi \\
z_1 &= z - \zeta \\
z_3 &= z + \zeta \\
\end{split}$$

$$\eq{
R = |x_1|, \quad v_1 = |z_1|, \quad v_3 = |z_3|
}$$

$$\begin{split}
d_1 &= R^2 \\
d_2 &= v_1^2 \\
d_3 &= v_3^2 \\
d_4 &= d_1 + d_2 \\
d_5 &= d_4^2 \\
d_6 &= d_1 + d_3 \\
d_7 &= d_6^2 \\
d_8 &= 2 R
\end{split}$$

$$\eq{
r_1 = \sqrt{d_4}, \quad r_3 = \sqrt{d_6}
}$$

$$\eq{
X = K R, \quad V_3 = K v_3, \quad Z = V_3 - \mathrm{i} X
}$$

$$\begin{split}
k_1 &= 2 K \\
k_2 &= k_1 K \\
\end{split}$$

$$\begin{split}
e_1 &= e^{-Z} \\
e_2 &= e_1 \mathrm{E}_1(-Z) \\
e_3 &= e_2 + \frac{1}{Z} \\
e_4 &= e_3 + \frac{1}{Z^2} \\
e_5 &= 2 \pi e_1 \\
e_6 &= K e_5 \\
e_7 &= K e_6 \\
\end{split}$$

$$\begin{split}
s_{x_1} &= \operatorname{sign}(x_1) \\
s_{z_1} &= \operatorname{sign}(z_1) \\
\end{split}$$

$$\eq{
G^{\text{deep}} = \log(K r_1) + \log(K r_3) - 2 [ \mathfrak{Re}(e_2) + \log(|Z|) ] - \mathrm{i} e_5
}$$

$$\eq{
G^{\text{deep}} = \log\left(\frac{r_1}{r_3}\right) - 2 \mathfrak{Re}(e_2) - \mathrm{i} e_5, \quad \text{for } X \gt 1
}$$

$$\eq{
G_x^{\text{deep}} = s_{x_1} \left( \frac{R}{d_4} - \frac{R}{d_6} + k_1 \mathfrak{Im}(e_3) + e_6 \right)
}$$

$$\eq{
G_z^{\text{deep}} = s_{z_1} \frac{v_1}{d_4} + \frac{v_3}{d_6} - k_1 \mathfrak{Re}(e_3) - \mathrm{i} e_6
}$$

$$\eq{
G_{xx}^{\text{deep}} =\frac{d_2 - d_1}{d_5} + \frac{d_1 - d_3}{d_7} + k_2 \mathfrak{Re}(e_4) + \mathrm{i} e_7
}$$

$$\eq{
G_{zz}^{\text{deep}} = -G_{xx}^{\text{deep}}
}$$

$$\eq{
G_{xz}^{\text{deep}} = -s_{x_1} \left( s_{z_1} v_1 \frac{d_8}{d_5} + v_3 \frac{d_8}{d_7} - k_2 \mathfrak{Im}(e_4) - e_7 \right)
}$$

## Plots
The expressions provided in the previous section were implemented in Julia and plotted with [Makie] for $K = 1$ and a source point located at $\mathrm{P} = (0, -2)$.

{{< figure src="images/G.png" alt="G" align="center" >}}

{{< figure src="images/Gx.png" alt="Gx" align="center" >}}

{{< figure src="images/Gz.png" alt="Gz" align="center" >}}

{{< figure src="images/Gxx.png" alt="Gxx" align="center" >}}

{{< figure src="images/Gzz.png" alt="Gzz" align="center" >}}

{{< figure src="images/Gxz.png" alt="Gxz" align="center" >}}

In future posts, this Green function will be applied for the solution of the radiation and diffraction problems of an oscillating floating body.

## References
1. Ed Mackay. 2021. The Green function for diffraction and radiation of regular waves by two-dimensional structures. European Journal of Mechanics - B/Fluids 87 (May 2021), 151–160. https://doi.org/10.1016/j.euromechflu.2021.01.012

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[Makie]: https://docs.makie.org/
[SymPy]: https://www.sympy.org
