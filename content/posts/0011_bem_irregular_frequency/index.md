---
date: '2025-08-23'
draft: false
title: 'Removing irregular frequencies of the linear wave-body problem'
author: 'Rodrigo Castro'
summary: 'Irregular frequencies of the linear wave-body problem are removed using the extended boundary condition method.'
tags: ['Boundary Element Method', 'Potential Flow']
---

## Introduction
In the [last post][bem_fsgfunc], the problem caused by irregular frequencies was observed on the solution of the radiation and diffraction problems of a floating body, using the Boundary Element Method (BEM) powered by the [free-surface Green function][fsgfunc]. In this post, the integral equations for the extended boundary condition method are presented and used to remove irregular frequencies from the solution of the wave-body problem. Additionally, an extra degree of freedom is considered and new results are displayed.

## Methods
First, the equations for the extended boundary condition method are displayed. Compared to the [previous post][bem_fsgfunc], not much changes on how the quantities of interest are calculated (radiation coefficients, exciting forces, etc.), so their expressions are not presented. As an extra, this chapter ends by describing the calculation of the wave reflection and transmission coefficients, something that was not evaluated previously.

### Extended boundary condition method
To remove irregular frequencies using the extended boundary condition method, a rigid free surface condition is imposed on the interior free surface $S_i$, depicted in the image below.

{{< figure src="images/cylinder.svg" alt="Cylinder" align="center" >}}

Adapting what was presented by Zhu (1994) to the two-dimensional problem, the integral equation becomes:

$$\eq{
\alpha \phi(\boldsymbol{\xi}) + 
\int_{S} \phi(\mathbf{x}) \frac{\partial G(\mathbf{x}, \boldsymbol{\xi})}{\partial \mathbf{n}(\mathbf{x})} \,dS(\mathbf{x}) = 
\int_{S} \frac{\partial \phi(\mathbf{x})}{\partial \mathbf{n}(\mathbf{x})} G(\mathbf{x}, \boldsymbol{\xi}) \,dS(\mathbf{x}),
}$$

where $S = S_b + S_i$, and

$$\eq{
\alpha = \begin{cases}
    -\pi & \text{if } \boldsymbol{\xi} \in S_b,\\
    +2 \pi & \text{if } \boldsymbol{\xi} \in S_i.
\end{cases}
}$$

### Reflection and transmission coefficients
To investigate how much energy is reflected or transmitted by the wave-body interaction, we use the following expressions for the reflection and transmission coefficients, respectively:

$$\eq{
K_r = \left\lvert \frac{\bar{\phi}^R + \bar{\phi}^D}{\bar{\phi}^I} \right\rvert _{x = -\infty}
}$$

$$\eq{
K_t = \left\lvert \frac{\bar{\phi}^R + \bar{\phi}^D + \bar{\phi}^I}{\bar{\phi}^I} \right\rvert _{x = +\infty}
}$$

## Results
The results for the wave-body interaction are presented below. Comparisons are made with the analytical and experimental data obtained by Vugts (1968) for the circular cylinder oscillating in sway and heave. It's possible to notice the absence of irregular frequencies and an overall good agreement between the BEM and the reference results.

{{< figure src="images/radiation_coefs_sway.svg" alt="Radiation coefficients sway" align="center" >}}

{{< figure src="images/radiation_coefs_heave.svg" alt="Radiation coefficients heave" align="center" >}}

{{< figure src="images/radiated_wave_amplitude.svg" alt="Radiated wave amplitude" align="center" >}}

{{< figure src="images/exciting_forces.svg" alt="Exciting forces" align="center" >}}

The next two figures shows the Response Amplitude Operator (RAO) for sway and heave and the animation for the cylinder-wave system at the frequency of maximum RAO, which is at $\omega^\ast \approx 0.86$.

{{< figure src="images/RAO.svg" alt="RAO's" align="center" >}}

{{< figure src="images/cylinder1.gif" alt="Animation 1" align="center" >}}

The following plot displays the reflection and transmission coefficients. Notice the interesting behavior for $\omega^\ast \approx 1.1$, where reflection is maximum and transmission is minimum. At this condition, almost all of the incident wave is reflected by the oscillating cylinder, as depicted by the animation that follows.

{{< figure src="images/wave_coefficients.svg" alt="Reflection and transmission coefficients" align="center" >}}

{{< figure src="images/cylinder2.gif" alt="Animation 2" align="center" >}}

## Conclusion
This post briefly presented the successful removal of irregular frequencies of the wave-body problem solved by the Boundary Element Method powered by the free-surface Green function. The code developed until now has demonstrated it has all it is needed to solve the most basic two-dimensional linear wave-body problems in waters of infinite-depth. Now it's time to study how this code can be modified to solve the finite-depth problem, which will increase its range of applications.

## References
1. Xuemei Zhu. 1994. Irregular frequency removal from the boundary integral equation for the wave body problem. Master's thesis. Massachusetts Institute of Technology, Cambridge, MA, USA. https://dspace.mit.edu/handle/1721.1/11691
2. J. H. Vugts. 1968. The hydrodynamic coefficients for swaying, heaving and rolling cylinders in a free surface. International Shipbuilding Progress, 15, 167 (1968), 251–276. https://doi.org/10.3233/ISP-1968-1516702

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[fsgfunc]: ../0009_2d_inf_depth_fsurface_gfunction/
[bem_fsgfunc]: ../0010_bem_fsurface_gfunction/
