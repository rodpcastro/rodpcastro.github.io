---
date: '2025-08-11'
draft: false
title: 'Solving the floating body radiation and diffraction problems with the BEM and the free-surface Green function'
author: 'Rodrigo Castro'
summary: 'Solution of the two-dimensional floating body radiation and diffraction problems with the Boundary Element Method and the free-surface Green function.'
tags: ['Boundary Element Method', 'Potential Flow']
---

## Introduction
In this article, the Green function described in the [previous post][gfunction] is put to work. Combined with the Boundary Element Method (BEM) implemented in [TwoDuBEM], the free-surface Green function will be used to solve the radiation and diffraction problems of a two-dimensional floating cylinder.

## Methods
The radiation problem was already described and solved in an [earlier post][cylinderpy]. Here, we add what concerns the diffraction problem, the floating body response and the new mesh.

### Diffraction problem
In the diffraction problem, we consider incident regular waves moving past the fixed body. In addition to satisfying Laplace's equation and the free-surface and radiation boundary conditions, the superposition of the known incident wave potential $\phi^I$ and the diffraction (or disturbance) potential $\phi^D$ must satisfy the boundary condition of zero normal velocity on the body, which can be written, after dropping the temporal term $e^{-\mathrm{i} \omega t}$, as

$$\eq{
\frac{\partial}{\partial n} (\bar{\phi}^I + \bar{\phi}^D) = 0, \quad \text{on} \quad S,
}$$

where $n$ is the normal vector to the body, pointing outside of the fluid domain, and $S$ is the submerged surface of the body. The complex amplitude of the incident wave potential for infinite-depth is defined by

$$\eq{
\bar{\phi}^I(x, z) = \frac{\mathrm{i} g}{\omega} e^{K z} e^{\mathrm{i} K x},
}$$

where $K = \omega^2 / g$. One of the main goals of solving the diffraction problem is obtaining the force exerted on the body. Similar to the radiation force, this is obtained by

$$\eq{
f_z = \mathrm{i} \rho \omega \int_{S} (\bar{\phi}^I + \bar{\phi}^D) n_z \, ds,
}$$

where $n_z$, for heave, is simply the $z$-component of the normal vector to the body.

For obtaining the force exerted on the body, it is interesting to note that the solution of the diffraction problem is not actually necessary. This can also be obtained from the solution of the radiation problem by taking advantage of the Haskind relation, as described by Newman (1962). For two-dimensional bodies with transverse symmetry, the amplitude of the force $f_z$ can be evaluated from the heave radiation damping as

$$\eq{
|f_z| = \zeta_a \sqrt{\frac{\rho g^2}{\omega} b_{z}},
}$$

where $\zeta_a$ is the incident wave amplitude, assumed 1 for simplicity. Expression $(4)$ will be used to validate the results obtained from the solution of the diffraction problem. 

In the results section, the force will be displayed in the following dimensionless form

$$
f^\ast_z = \frac{|f_z|}{2 \rho g R \zeta_a}
$$

### Boundary Element Method
In the Boundary Element Method, the system of equations that solves for the diffraction potential $\bar{\phi}^D$ comes from equation $(1)$ and is simply

$$\eq{
\left[ \mathbb{Q} \right] \{ \bar{\phi}^D \} = -\left[ \mathbb{G} \right] \{ \bar{\phi}^I_n \},
}$$

where $\bar{\phi}^I_n$ is the derivative of $\phi^I$ in the direction of the normal to body's surface.

### Response Amplitude Operator
The floating body response to incident regular waves can be calculated for heave, after solving for the radiation and diffraction problems, as

$$\eq{
\bar{z} = \frac{f_z}{-\omega^2 (M + a_z) - \mathrm{i} \omega b_z + c_z},
}$$

where $M$ is mass of the body, $a_z$ and $b_z$ are the radiation coefficients, and $c_z$ is the restoring coefficient, which in two dimensions is simply

$$\eq{
c_z = \rho g l_{w},
}$$

and $l_w$ is the body's waterline length, which for the cylinder is $2 R$.

The Response Amplitude Operator (RAO) and the phase angle in radians for heave can be computed by

$$\eq{
RAO_{heave} = |\bar{z}|, \quad \epsilon = \arctan{\frac{\mathfrak{Im}(\bar{z})}{\mathfrak{Re}(\bar{z})}}
}$$

### Mesh
Since the infinite-depth free-surface Green function already satisfies the boundary conditions at the free surface and the radiation condition at infinity, all that remains is the boundary condition on the body's surface. This means that the boundary mesh is reduced to what is depicted in the following image, where 10 elements are enough to represent the floating cylinder.

{{< figure src="images/mesh.svg" alt="Mesh" align="center" >}}

## Results
First, the results obtained from the solution of the radiation problem are presented again, but this time, the BEM was powered by the free-surface Green function. The analytical (Ursell) and experimental (Vugts) results used as references were discussed in a [previous article][cylinderjl]

All results obtained with the BEM contain a singular behavior for $\omega^\ast \approx 1.36$. Frequencies where this happens are known in the literature as [irregular frequencies][irrfreq], which are not associated with a physical wave problem, but are a mathematical and computational technicality. This problem will be further discussed and tackled in a future post.

Ignoring the errors introduced by the irregular frequency phenomena, the solution of the radiation problem provides accurate results when compared with both analytical and experimental results.

{{< figure src="images/heave_added_mass.svg" alt="Added mass" align="center" >}}

{{< figure src="images/heave_wave_damping.svg" alt="Wave damping" align="center" >}}

{{< figure src="images/wave_rao.svg" alt="Wave RAO" align="center" >}}

Ursell did not provide expressions for the exciting force in heave, but this can be computed through Haskind relation, as explained [previously](#diffraction-problem). The next figure shows how the force in heave obtained through BEM also matches analytical and experimental data.

{{< figure src="images/heave_force.svg" alt="Force" align="center" >}}

The response of the floating cylinder is computed according to $(5)$ and presented below for reference. The RAO approaching 1 for small frequencies is expected, since this is the condition where the floating cylinder heaves in phase with the wave and with the same amplitude.

{{< figure src="images/heave_RAO.svg" alt="Heave RAO" align="center" >}}

As a cherry on the top, we end this results section with two animations of the heaving cylinder, the final result after solving the radiation and diffraction problems successfully. The first animation is for the frequency when the RAO is maximum.

{{< figure src="images/cylinder1.gif" alt="Cylinder animation" align="center" >}}

The second animation is for $\omega^\ast \approx 1.22$ (**not** in the region strongly affected by the irregular frequency). At this frequency, it's possible to see that the wave amplitude on the left side of the cylinder is higher than the wave amplitude on the right side, which means that the cylinder is acting as a barrier. If the horizontal movement had also been considered, this difference could have been less pronounced, but this is a task for another article.

{{< figure src="images/cylinder2.gif" alt="Cylinder animation" align="center" >}}

## Conclusion
This post showed once again the sucessful use of the BEM for solving potential flow problems. This floating cylinder might appear one more time when the removal of irregular frequencies is achieved.

## References
1. J. N. Newman. 1962. The exciting forces on fixed bodies in waves. J. Ship Res. 6, 4 (1962), 10–17. https://doi.org/10.5957/jsr.1962.6.4.10

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[cylinderjl]: ../0006_julia_2d_cylinder_free_surface/
[cylinderpy]: ../0008_bem_floating_cylinder/
[gfunction]: ../0009_2d_inf_depth_fsurface_gfunction/
[twodubem]: ../../projects/twodubem/
[irrfreq]: https://www.orcina.com/webhelp/OrcaWave/Content/html/Theory,Irregularfrequencies.htm
