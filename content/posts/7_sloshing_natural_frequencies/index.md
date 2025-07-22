---
date: '2025-07-22'
draft: false
title: 'Natural frequencies and modes of sloshing in a rectangular tank'
author: 'Rodrigo Castro'
summary: 'The natural frequencies and modes of sloshing in a rectangular tank are obtained through solution of a potential flow.'
tags: ['Sloshing', 'Boundary Element Method', 'Python', 'Potential Flow']
---

## Introduction
This post displays the computation of natural frequencies and modes of sloshing in a rectangular tank by using the two-dimensional boundary element method implemented in [TwoDuBEM]. In addition to the numerical method, the analytical solution is briefly presented. Experimental data is also displayed in the results section.

## Methods
The problem to be solved in this post is the sloshing in a rectangular container. The following image depicts the tank and the relevant dimensions of the problem: the tank's width (or breadth) $L$ and its depth $h$. The coordinate system is defined such that its origin is placed at the midpoint between the sidewalls, the $x$ axis is coincident with the undisturbed free surface while the $y$ axis is normal to it. The fluid domain is enclosed by two boundaries: the free surface $S_F$ and the rigid boundary of the tank $S_R$. Unit normal vectors to the boundaries point outward from the fluid domain.

{{< figure src="images/tank.svg" alt="Tank" align="center" >}}

Assuming that the flow inside the tank is incompressible and irrotational, the equation governing the flow is the Laplace equation:

$$\eq{
\nabla^2 \phi = 0,
}$$

where $\phi$ is the velocity potential.

If the container is stationary while the liquid is free to oscilate in simple harmonic motion with frequency $\omega$, then all time-dependent quantities may me assumed to be varying simple harmonically as well. $\phi$ is expressed as

$$\eq{
\phi(x, y, t) = \bar{\phi}(x, y) e^{i \omega t},
}$$

and as consequence, $(1)$ can be rewritten to

$$\eq{
\nabla^2 \bar{\phi} = 0.
}$$

The boundary condition on the walls of the tank $S_R$ is 

$$\eq{
\frac{\partial \bar{\phi}}{\partial n} = 0,
}$$

and on the liquid free surface $S_F$, the linearized condition is

$$\eq{
\frac{\partial \bar{\phi}}{\partial n} = \frac{\omega^2}{g} \bar{\phi},
}$$

where $g$ is the acceleration of gravity.

In the following section, the numerical approach of the boundary element method implemented in TwoDuBEM is revisited and extended to solve the present problem.

### Numerical
On previous posts on the boundary element method ([post 1][bem_python], [post 2][bem_element]), it was shown how a differential equation could be transformed into a boundary integral equation, and how the latter could be approximated by a system of equations of the form:

$$\eq{
\frac{1}{2} \bar{\phi}(\mathbf{x}_i) = 
\sum_{j=1}^{N} \bar{\phi}_j \mathbb{Q}(\mathbf{x}_j, \mathbf{x}_i) -
\sum_{j=1}^{N} \bar{q}_j \mathbb{G}(\mathbf{x}_j, \mathbf{x}_i),
}$$

where

$$\eq{
\bar{q}_j = \frac{\partial \bar{\phi}_j}{\partial n},
}$$

$$\eq{
\mathbb{G}(\mathbf{x}_j, \mathbf{x}_i) = 
\int_{S_j} \frac{1}{2 \pi} \ln |\mathbf{x} - \mathbf{x}_i| \,dS,
}$$

$$\eq{
\mathbb{Q}(\mathbf{x}_j, \mathbf{x}_i) = 
\int_{S_j} \frac{\partial}{\partial n} \left( \frac{1}{2 \pi} \ln |\mathbf{x} - \mathbf{x}_i| \right) \,dS.
}$$

In matrix notation, $(6)$ can be written as

$$\eq{
\left[ \mathbb{Q} \right] \{ \bar{\phi} \} = \left[ \mathbb{G} \right] \{ \bar{q} \},
}$$

where

$$\eq{
\mathbb{Q}_{ij} = \begin{cases}
    \mathbb{Q}(\mathbf{x}_j, \mathbf{x}_i) & \text{if } i \neq j,\\
    \mathbb{Q}(\mathbf{x}_j, \mathbf{x}_i) - \frac{1}{2} & \text{if } i = j.
\end{cases}
}$$

Following the procedure described by Dutta (2000), equation $(10)$ can be arranged in the following form

$$\eq{
\left[ \mathbb{Q}_R \,|\, \mathbb{Q}_F \right]
\begin{Bmatrix}
    \bar{\phi}_R \\
    \bar{\phi}_F \\
\end{Bmatrix} = 
\left[ \mathbb{G}_R \,|\, \mathbb{G}_F \right]
\begin{Bmatrix}
    \bar{q}_R \\
    \bar{q}_F \\
\end{Bmatrix},
}$$

where the subscripts $R$ and $F$ denote the rigid boundary and the free surface, respectively. In this matrix representation, the values corresponding to the free surface are arranged in the right-most part of $[\mathbb{Q}]$ and $[\mathbb{G}]$ and on the bottom part of each column vector $\{\bar{\phi}\}$ and $\{\bar{q}\}$.

Now, we do the following rearrangement

$$\eq{
\left[ \mathbb{Q}_R \,|\, -\mathbb{G}_F \right]
\begin{Bmatrix}
    \bar{\phi}_R \\
    \bar{q}_F \\
\end{Bmatrix} = 
\left[ \mathbb{G}_R \,|\, -\mathbb{Q}_F \right]
\begin{Bmatrix}
    \bar{q}_R \\
    \bar{\phi}_F \\
\end{Bmatrix},
}$$

and substituting boundary conditions $(4)$ and $(5)$, we end up with

$$\eq{
\left[ \mathbb{Q}_R \,|\, -\mathbb{G}_F \right]
\begin{Bmatrix}
    \bar{\phi}_R \\
    \bar{q}_F \\
\end{Bmatrix} = 
\left[ \mathbb{G}_R \,|\, -\mathbb{Q}_F \right]
\begin{Bmatrix}
    0 \\
    \frac{g}{\omega^2} \bar{q}_F \\
\end{Bmatrix}.
}$$

Now, we define the following matrix

$$\eq{
[\mathbb{C}] = 
\begin{bmatrix}
    0 & 0 \\
    0 & I_F \\
\end{bmatrix},
}$$

where $I_F$ is the identity matrix of size $N_F \times N_F$, and $N_F$ is the number of elements used on the free surface $S_F$. Using the matrix $\mathbb{C}$, we rewrite $(14)$ as

$$\eq{
\left[ \mathbb{Q}_R \,|\, -\mathbb{G}_F \right]
\begin{Bmatrix}
    \bar{\phi}_R \\
    \bar{q}_F \\
\end{Bmatrix} = 
\left[ \mathbb{G}_R \,|\, -\mathbb{Q}_F \right]
[\mathbb{C}]
\begin{Bmatrix}
    \bar{\phi}_R \\
    \bar{q}_F \\
\end{Bmatrix}
\frac{g}{\omega^2}.
}$$

Defining

$$\eq{
[\mathbb{T}] = [\mathbb{Q}_R \,|\, -\mathbb{G}_F]^{-1} [\mathbb{G}_R \,|\, -\mathbb{Q}_F] [\mathbb{C}],
}$$

we can rewrite equation $(16)$ as the following eigenvalue problem

$$\eq{
[ \mathbb{T}]
\begin{Bmatrix}
    \bar{\phi}_R \\
    \bar{q}_F \\
\end{Bmatrix} = 
\frac{\omega^2}{g}
\begin{Bmatrix}
    \bar{\phi}_R \\
    \bar{q}_F \\
\end{Bmatrix}.
}$$

Knowing that the first columns of $\mathbb{T}$ corresponding to the rigid boundary are all null, equation $(18)$ can be further simplified to

$$\eq{
[T_{FF}]\{ \bar{q}_F \} = \frac{\omega^2}{g} \{ \bar{q}_F \},
}$$

from which the natural frequencies of oscillation and the natural sloshing modes of the free surface may be found.

The prodecure described in this section is implemented as an extension of TwoDuBEM, without changing the original code. See the files in the [appendix](#appendices) for more details.

### Analytical
The analytical expressions for the natural frequencies and modes can be found on Ibrahim (2005). First, defining the wave number $k_m$ for the mode $m$ as 

$$\eq{
k_m = \begin{cases}
\frac{2 m \pi}{L} & \text{for symmetric modes}, \\
\frac{(2 m - 1) \pi}{L} & \text{for asymmetric modes},
\end{cases}
}$$

the natural frequencies of the free surface are given by the expression

$$\eq{
\omega_m^2 = g k_m \tanh(k_m h).
}$$

Symmetric mode shapes have the form of $\cos(k_m x)$ and asymmetric mode shapes have the form of $\sin(k_m x)$.

## Results
The first set of results to be presented are the natural modes. The numerical output was obtained for container of dimensions $L=1.0$ and $h=0.5$, and a mesh of 23 elements (5 along the bottom, 4 along the sides and 10 on the free surface). It's noticeable that even with a small number of elements on the free surface, the numerical mode shapes follow the analytical reference very precisely. Of course, without the analytical reference, these modes would have to be obtained with many more elements on the free surface.

{{< figure src="images/natural_modes.svg" alt="Natural modes" align="center" >}}

Next figure displays the analytical and numerical natural frequencies for the first 4 natural modes. Numerical results were derived with a mesh of 46 elements (10 along the bottom, 8 along the sides and 20 on the free surface). Notice how the difference increases as the frequency gets higher. This makes sense, because higher frequency means larger variation. In order to capture this variation, the mesh needs to be further refined.

{{< figure src="images/natural_frequencies.svg" alt="Natural frequencies" align="center" >}}

Last figure is a special case of the last result, but with the addition of the experimental data extracted from the report published by Addington (1960). From this figure, we can be certain that the mathematical model described earlier is a good representation of the real phenomena. A larger difference is observed between theoretical and experimental results in the realm of higher natural frequencies, indicating a region where fluid viscosity begins to play a role.

{{< figure src="images/natural_frequency_mode1.svg" alt="Natural frequency mode 1" align="center" >}}

## Conclusion
This was another sucessfull application of the boundary element method for the solution of potential flow. The next application on the topic of sloshing will be the evaluation of the forces due to lateral motion (surge) or rotation (pitch) of the tank. 

## References
1. S. Dutta and M. K. Laha. 2000. Analysis of the small amplitude sloshing of a liquid in a rigid container of arbitrary shape using a low-order boundary element method. Int. J. Numer. Meth. Engng. 47, 9 (March 2000), 1633–1648. https://doi.org/10.1002/%28SICI%291097-0207%2820000330%2947%3A9%3C1633%3A%3AAID-NME851%3E3.0.CO%3B2-1
2. R. A. Ibrahim. 2005. Liquid Sloshing Dynamics: Theory and Applications. Cambridge University Press, Cambridge.
3. J. W. Addington. 1960. Dynamics of Fuel in Tanks. Report CoA/N-99. College of Aeronautics, Cranfield, Bedfordshire. https://reports.aerade.cranfield.ac.uk/handle/1826.2/4460

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[bem_python]: ../2_bem_python/
[bem_element]: ../3_2d_constant_boundary_element/
[twodubem]: https://github.com/rodpcastro/twodubem
