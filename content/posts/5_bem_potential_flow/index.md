---
date: '2025-07-11'
draft: false
title: 'Solving the 2D potential flow using the Boundary Element Method'
author: 'Rodrigo Castro'
summary: 'Computation of the 2D velocity potential around a circular cylinder using the Boundary Element Method.'
tags: ['Boundary Element Method', 'Potential Flow']
---

## Introduction
In the conclusion section of the [first post about the boundary element method][bem_python], one of the recommended improvements is the refactoring of the code using a object oriented approach. The past two weeks have been dedicated to this task through the creation of the [TwoDuBEM] project.

In this post, the goal is to show the use of the [Boundary Element Method][bem] implemented in TwoDuBEM to solve the two-dimensional [potential flow][potflow] around a circular cylinder placed in a uniform flow.

## Methods
The potential flow is described by the velocity potential that satisfies the Laplace's equation, which fundamental solution is already known and studied in a [previous post][2d_incoef]. The following topics describe the boundary value problem that represents the potential flow around a circular cylinder and the analytical solution to that problem.

## Boundary value problem
The velocity potential $\Phi$ describing the potential flow around a circular cylinder of radius $R$ placed in a uniform flow of speed $U$ satisfies:

$$\eq{
\nabla^2 \Phi = 0 \quad \text{for} \quad r = \sqrt{x^2 + y^2} \ge R,
}$$

where $\Phi$ is the composition of the potential due to a uniform flow and the unkown distubance potential $\phi$.

$$\eq{
\Phi = Ux + \phi.
}$$

The no-penetration boundary condition states that:

$$\eq{
\frac{\partial \Phi}{\partial n} = 0 \quad \text{for} \quad \sqrt{x^2 + y^2} = R,
}$$

where $n$ is the unit vector normal to the boundary, which points outward from the domain in TwoDuBEM's reference. In the case of a circular cylinder, this vector is coincident with the radial direction, but with opposite sign, and since $x$ can be written in polar coordinates as

$$ x = r \cos\theta, $$

where $r = \sqrt{x^2 + y^2}$ and $\theta$ is defined as the angle between $r$ and the positive $x$ axis. Therefore, the left side of $(3)$ can be written as:

$$\eq{
\frac{\partial \Phi}{\partial n} = -U \cos\theta + \frac{\partial \phi}{\partial n},
}$$

and the boundary condition for the unkown potential $\phi$ is:

$$\eq{
\frac{\partial \phi}{\partial n} = U \frac{x}{R} \quad \text{for} \quad \sqrt{x^2 + y^2} = R.
}$$

## Analytical solution

The analytical solution for the total potential $\Phi$ is equivalent to the composition of the uniform flow and dipole potentials. The dipole intensity can be calculated as a function of the uniform flow speed and the cylinder radius, giving the following expressions for the total potential and its gradient in cartesian coordinates:

$$\eq{
\Phi = U x \left( 1 + \frac{R^2}{x^2 + y^2} \right),
}$$

$$\eq{
\Phi_x = U \left( 1 + \frac{R^2}{x^2 + y^2} - \frac{2 R^2 x^2}{(x^2 + y^2)^2} \right),
}$$

$$\eq{
\Phi_y = \frac{-2 U R^2 x y}{(x^2 + y^2)^2}.
}$$

The tangential velocity at the boundary of the cylinder as a function of the angle $\theta$ is given by:

$$\eq{
V = 2 U \sin\theta.
}$$

From Bernoulli's theorem $(10)$ and the tangential velocity in $(9)$, we can derive the pressure coefficient $C_p$:

$$ p_0 - \frac{1}{2} \rho U^2 = p - \frac{1}{2} \rho V^2, \tag{10} $$

$$ C_p = \frac{p - p_0}{\frac{1}{2} \rho U^2} = 1 - 4 \sin^2\theta. \tag{11} $$


## Results
In the following set of images, the numerical results were obtained with 360 elements. The large number of elements was necessary to achieve better accuracy on the computation of the gradient.

The first image displays the results for the total potential $\Phi$ and its gradient components in cartesian coordinates.

{{< figure src="images/total_potential.svg" alt="Total potential" align="center" >}}

The next image shows the streamlines around the cylinder.

{{< figure src="images/streamlines.svg" alt="Total potential" align="center" >}}

This last image shows the anaytical and numerical results for the pressure coefficient $C_p$.

{{< figure src="images/pressure_coefficient.svg" alt="Total potential" align="center" >}}

## Conclusion
This post was a demonstration of [TwoDuBEM] in action. Future posts might bring significant improvements to the project or other applications of the implemented algorithm.

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[bem_python]: ../2_bem_python/
[2d_incoef]: ../3_2d_constant_boundary_element/
[twodubem]: https://github.com/rodpcastro/twodubem
[bem]: https://en.wikipedia.org/wiki/Boundary_element_method 
[potflow]: https://en.wikipedia.org/wiki/Potential_flow
