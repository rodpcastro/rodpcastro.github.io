---
date: '2025-05-17'
draft: true
title: '2D constant boundary element influence coefficients'
author: 'Rodrigo Castro'
summary: 'Computation of influence coefficients for 2D constant boundary elements.'
tags: ['Boundary Element Method']
---

## Introduction
This post presents a first improvement on the [previous post][2_bem_python], specially on the calculation of the influence coefficients for the two dimensional constant boundary element.

The ojective is to study how the [scipy.integrate.quad][scipy_quad] function and the [Gauss-Legendre] quadrature compare to the analytical evaluation of the integrals that represent the influence coefficients.

## Methods
The [previous post][2_bem_python] introduced the boundary integral equation for the Laplace's equation. The equation is presented below with different symbols to match the nomenclature used by *Crouch* and *Mogilevskaya* in their [book][fcbem], which will be the main reference for this post.

$$\eq{
\frac{1}{2} \phi(\boldsymbol{\xi}) = 
\int_{S} \phi(\mathbf{x}) Q(\mathbf{x}, \boldsymbol{\xi}) \,dS(\mathbf{x}) -
\int_{S} q(\mathbf{x}) G(\mathbf{x}, \boldsymbol{\xi}) \,dS(\mathbf{x}).
}$$

The fundamental solution and their normal derivative are:

$$\eq{
G(\mathbf{x}, \boldsymbol{\xi}) = \frac{1}{2 \pi} \ln |\mathbf{x} - \boldsymbol{\xi}|,
}$$

$$\eq{
Q(\mathbf{x}, \boldsymbol{\xi}) = \frac{\partial G(\mathbf{x}, \boldsymbol{\xi})}{\partial \mathbf{n}}.
}$$

After discretizing the boundary $S$ into $N$ straight elements with constant $\phi_j$ and $q_j$, the boundary integral equation in $(1)$ can be approximated by

$$\eq{
\frac{1}{2} \phi(\boldsymbol{\xi}) = 
\sum_{j=1}^{N} \phi_j \mathbb{Q}(\mathbf{x}_j, \boldsymbol{\xi}) -
\sum_{j=1}^{N} q_j \mathbb{G}(\mathbf{x}_j, \boldsymbol{\xi}),
}$$

where

$$\eq{
\mathbb{G}(\mathbf{x}_j, \boldsymbol{\xi}) = 
\int_{S_j} G(\mathbf{x}, \boldsymbol{\xi}) \,dS(\mathbf{x}),
}$$

$$\eq{
\mathbb{Q}(\mathbf{x}_j, \boldsymbol{\xi}) = 
\int_{S_j} Q(\mathbf{x}, \boldsymbol{\xi}) \,dS(\mathbf{x}).
}$$

The following subtopics present how to evaluate the integrals $(5)$ and $(6)$ by three methods. The first method is analytical and the other two are numerical, using the Gauss-Legendre quadrature and the `scipy.integrate.quad` function.

### Analytical integration

{{< figure src="images/intGQ.svg" alt="analytical G Q parameters" align="center">}}

### Gauss-Legendre quadrature

### Scipy.integrate.quad

## Results


## Conclusion


## References
1. Steven L. Crouch, Sofia G. Mogilevskaya. 2024. A First Course in Boundary Element Methods. Springer, Cham, Switzerland.

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[2_bem_python]: ../2_bem_python
[scipy_quad]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
[Gauss-Legendre]: https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
[fcbem]: https://search.worldcat.org/title/1450559181
