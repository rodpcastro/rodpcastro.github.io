---
date: '2025-08-01'
draft: true
title: 'Using the Boundary Element Method to find the radiation potential of a floating cylinder'
author: 'Rodrigo Castro'
summary: 'Radiation coefficients and the wave amplitude are computed for a cylinder heaving on the free surface of a liquid.'
tags: ['Boundary Element Method', 'Python', 'Potential Flow']
---

## Introduction
It is time to revisit a problem solved in an [earlier post][j2d_cylinder]: a cylinder oscilating in heave on the free surface of a liquid. The difference now is the solution method used, which is based on the Boundary Element Method implemented in [TwoDuBEM].

## Methods
The next subtopics describe the problem and the discretization for the boundary element method.

### Formulation of the problem
The following image depicts the cylinder and boundaries of the fluid domain.

{{< figure src="images/fluid_domain.svg" alt="Fluid domain" align="center">}}

Assuming harmonic motion of the cylinder with frequency $\omega$, all time time-dependent quantities can be assumed to be varying harmonically with the same frequency. Then, the radiation potential can be expressed as

$$\eq{
\phi = \bar{\phi} e^{\mathrm{i} \omega t},
}$$

and the problem can be solved in the frequency domain in terms of the complex amplitude of the radiation potential $\bar{\phi}$, which must satisfy Laplace's equation in the fluid domain

$$\eq{
\nabla^2 \bar{\phi} = 0, 
}$$

the boundary condition on the bottom

$$\eq{
\frac{\partial \bar{\phi}}{\partial z} = 0 \quad \text{on} \quad S_B, 
}$$

the free surface boundary condition

$$\eq{
\frac{\partial \bar{\phi}}{\partial z} = \frac{\omega^2}{g} \bar{\phi} \quad \text{on} \quad S_F, 
}$$

the [Sommerfeld radiation condition][sommerfeld]

$$\eq{
\frac{\partial \bar{\phi}}{\partial \lvert x \rvert} + \mathrm{i} k \bar{\phi} \to 0, \quad \text{as} \quad \lvert x \rvert \to \infty, }$$

where the wavenumber $k$ is related to the wave frequency $\omega$ and the water depth $h$ by the dispersion relation

$$\eq{
\omega^2 = k g \tanh(k h)
}$$

Finally, the boundary condition on the surface of the cylinder is given by

$$\eq{
\frac{\partial \bar{\phi}}{\partial n} = \mathrm{i} \omega h_z \quad \text{on} \quad S_C, 
}$$

where $h_z$ is the heave mode shape.

After solving for $\bar{\phi}$, the hydrodynamic force in heave can be computed as

$$\eq{
F_z = \mathrm{i} \rho \omega \int_{S_C} \bar{\phi} n_z \, ds,
}$$

where $n_z$ is the component of the normal vector on the surface of the cylinder $S_C$ in the direction of the mode shape $h_z$. Then, the added mass $a_z$ and wave damping $b_z$ in heave are calculated as

$$\eq{
a_z = -\frac{\mathfrak{Re}(F_z)}{\omega^2}
}$$

$$\eq{
b_z = \frac{\mathfrak{Im}(F_z)}{\omega}
}$$

The complex wave amplitude $\bar{\zeta}$ at a distance from the floating cylinder is constant and given by

$$\eq{
\bar{\zeta} = -\frac{\mathrm{i} \omega}{g} \bar{\phi}
}$$

### Boundary Element Method
In the Boundary Element Method the system of equations to be solved, in matrix form, is

$$\eq{
\left[ \mathbb{Q} \right] \{ \bar{\phi} \} = \left[ \mathbb{G} \right] \{ \bar{q} \},
}$$

where

$$\eq{
\bar{q} = \frac{\partial \bar{\phi}}{\partial n},
}$$

and $\mathbb{G}$ and $\mathbb{Q}$ are the [influence matrices][2d_cbem].

Equation $(12)$ can be expanded to show submatrices corresponding to each boundary

$$\eq{
\left[ \mathbb{Q}_F \,|\, \mathbb{Q}_R \,|\, \mathbb{Q}_B \,|\, \mathbb{Q}_C \right]
\begin{Bmatrix}
    \bar{\phi}_F \\
    \bar{\phi}_R \\
    \bar{\phi}_B \\
    \bar{\phi}_C \\
\end{Bmatrix} = 
\left[ \mathbb{G}_F \,|\, \mathbb{G}_R \,|\, \mathbb{G}_B \,|\, \mathbb{G}_C \right]
\begin{Bmatrix}
    \bar{q}_F \\
    \bar{q}_R \\
    \bar{q}_B \\
    \bar{q}_C \\
\end{Bmatrix},
}$$

where the subscripts indicate a boundary: $F$ - free surface, $R$ - radiation surface, $B$ - bottom, $C$ - cylinder. Each of the submatrices of $\mathbb{G}$ and $\mathbb{Q}$ have $N$ lines and the number of columns equal the number of panels used to represent the corresponding boundary, and $N$ is total number of elements used.

Substituting the boundary conditions described in the [prior subtopic](#formulation-of-the-problem) in equation $(14)$, we get

$$\eq{
\left[ \mathbb{Q}_F \,|\, \mathbb{Q}_R \,|\, \mathbb{Q}_B \,|\, \mathbb{Q}_C \right]
\begin{Bmatrix}
    \bar{\phi}_F \\
    \bar{\phi}_R \\
    \bar{\phi}_B \\
    \bar{\phi}_C \\
\end{Bmatrix} = 
\left[ \mathbb{G}_F \,|\, \mathbb{G}_R \,|\, \mathbb{G}_B \,|\, \mathbb{G}_C \right]
\begin{Bmatrix}
    \frac{\omega^2}{g} \bar{\phi}_F \\
    -\mathrm{i} k \bar{\phi}_R \\
    0 \\
    \mathrm{i} \omega h_z \\
\end{Bmatrix}.
}$$

Rearranging the expression above, we get the system of equations that will give us the amplitude of the radiation potential $\bar{\phi}$ on each boundary

$$\eq{
\left[ \mathbb{Q}_F - \frac{\omega^2}{g} \mathbb{G}_F \,|\, \mathbb{Q}_R + \mathrm{i} k \mathbb{G}_R \,|\, \mathbb{Q}_B \,|\, \mathbb{Q}_C \right]
\begin{Bmatrix}
    \bar{\phi}_F \\
    \bar{\phi}_R \\
    \bar{\phi}_B \\
    \bar{\phi}_C \\
\end{Bmatrix} =  \\
\left[ \mathbb{G} \right]
\begin{Bmatrix}
    0 \\
    0 \\
    0 \\
    \mathrm{i} \omega h_z \\
\end{Bmatrix}.
}$$

The hydrodynamic force $F_z$ in equation $(8)$, from which the added mass and wave damping are computed, is calculated on the cylinder surface as

$$\eq{
F_z = \mathrm{i} \rho \omega \sum_{i=1}^{N_C} \bar{\phi}^i_C n^i_z s_i.
}$$

The wave amplitude is calculated as the mean value of expression $(11)$ for collocation points on the free surface situated at a distance larger than or equal to $\lambda + R$, where $\lambda$ is the wavelength and $R$ is the cylinder radius.

### Boundary mesh
The mesh for the Boundary Element Method was built as a function of the wavelength. This was done, because the range of wavelengths studied is too large, and a single mesh to handle the full range is not recommended as the number of elements can get unecessarily large. Alternatively, the mesh could be constructed for subranges, but experiments suggested that rebuilding the mesh, and the influence matrices, for each wavelength studied did not impose too much of a computational burden.
The mesh is constructed in such a way that the element size on the surface of the cylinder is $\lambda/20$ and grows as a geometric prossion the further it gets away from the cylinder, in the $x$ and $z$ directions. The image below presents the mesh used for one of the wavelengths studied.

{{< figure src="images/mesh.svg" alt="Mesh" align="center" >}}

## Results

{{< figure src="images/heave_added_mass.svg" alt="Added mass" align="center" >}}

{{< figure src="images/heave_wave_damping.svg" alt="Wave damping" align="center" >}}

{{< figure src="images/wave_rao.svg" alt="Wave RAO" align="center" >}}

## Conclusion


## References
1. Gabriel Maarten and Peter Wellens. 2022. A two-dimensional boundary element method with generating absorbing boundary condition for floating bodies of arbitrary shape in the frequency domain. International Shipbuilding Progress 69, 2 (Feb. 2022), 139–159. https://doi.org/10.3233/ISP-210007.

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[j2d_cylinder]: ../6_julia_2d_cylinder_free_surface/
[twodubem]: ../../projects/twodubem/
[sommerfeld]: https://www.wikiwaves.org/index.php/Sommerfeld_Radiation_Condition
[2d_cbem]: ../3_2d_constant_boundary_element/
