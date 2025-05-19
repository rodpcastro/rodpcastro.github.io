---
date: '2025-05-17'
draft: true
title: '2D constant boundary element influence coefficients'
author: 'Rodrigo Castro'
summary: 'Computation of influence coefficients for 2D constant boundary elements.'
tags: ['Boundary Element Method']
---

## Introduction
This post presents a first improvement on the [previous article][2_bem_python], specially on the calculation of the influence coefficients for the two dimensional boundary element method.

The ojective is to study how the [`scipy.integrate.quad`][spquad] function and the [Gauss-Legendre quadrature][GLquad] compare to the analytical evaluation of the integrals that represent the influence coefficients.

## Methods
The [previous post][2_bem_python] introduced the boundary integral equation for the Laplace's equation. The equation is presented below with a similar nomenclature to that used by *Crouch* and *Mogilevskaya* in their [book][fcbem], which will be the main reference for this work.

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
Integrals $(5)$ and $(6)$ are evaluated over a straight boundary element $S_j$ of length $2a_j$, shown as a red color line in the following figure. The element local system has the $x$ axis colinear with the element and the $y$ axis normal to it. The blue dot represents the source point $\boldsymbol{\xi}$, which has coordinates $(\xi_x, \xi_y)$ in the element local system. The other displayed paramaters will help with the analytical integration and are given in equations $(7)$ and $(8)$.

{{< figure src="images/intGQ.svg" alt="analytical G Q parameters" align="center">}}

$$\eq{
r_1 = \sqrt{\left(x - a_j \right)^2 + \xi_y^2}, \quad r_2 = \sqrt{\left(x + a_j \right)^2 + \xi_y^2},
}$$

$$\eq{
\theta_1 = \arctan \frac{\xi_y}{\xi_x - a_j}, \quad \theta_2 = \arctan \frac{\xi_y}{\xi_x + a_j}.
}$$

$\theta_1$ and $\theta_2$ can assume values between $-\pi$ and $\pi$.

With these new definitions, integrals $(5)$ and $(6)$ are rewritten as

$$\eq{
\mathbb{G}(\mathbf{x}_j, \boldsymbol{\xi}) = 
\int_{-a_j}^{a_j} G(x - \xi_x, \xi_y) \,dx,
}$$

$$\eq{
\mathbb{Q}(\mathbf{x}_j, \boldsymbol{\xi}) = 
\int_{-a_j}^{a_j} Q(x - \xi_x, \xi_y) \,dx,
}$$

where $G$ and $Q$ are

$$\eq{
G(x - \xi_x, y - \xi_y) = 
\frac{1}{2 \pi} \ln \sqrt{\left( x - \xi_x \right)^2 + \left( y - \xi_y \right)^2},
}$$

$$\eq{
Q(x - \xi_x, y - \xi_y) = 
\frac{1}{2 \pi} \frac{y - \xi_y}{\left( x - \xi_x \right)^2 + \left( y - \xi_y \right)^2}.
}$$

Expressions $(9)$ and $(10)$ can be analytically evaluated to give

$$\eq{
\mathbb{G}(\mathbf{x}_j, \boldsymbol{\xi}) = 
\frac{1}{2 \pi} \left[ \xi_y (\theta_1 - \theta_2) - (\xi_x - a_j) \ln r_1 + (\xi_x + a_j) \ln r_2 - 2a_j \right],
}$$

$$\eq{
\mathbb{Q}(\mathbf{x}_j, \boldsymbol{\xi}) = 
\frac{1}{2 \pi} (\theta_1 - \theta_2).
}$$

For any $\boldsymbol{\xi}$ colinear with the element, $\mathbb{Q}$ is evaluated to

$$\eq{
\mathbb{Q}(\mathbf{x}_j, \boldsymbol{\xi}) = 0.
}$$

The value of $\mathbb{G}$ for $\boldsymbol{\xi}$ coincident with the element's endpoints $(\pm a_j, 0)$ is obtained considering the limit case of the expression $(13)$ as $\xi_y = 0$ and $\xi_x \to \pm a_j$. 

$$\eq{
\mathbb{G}(\mathbf{x}_j, (\pm a_j, 0)) = \frac{a_j}{\pi} [\ln(2 a_j) - 1].
}$$

Similarly, the diagonal (self-effect) terms in $(4)$ can be obtained from the limit case of $(13)$ as $\boldsymbol{\xi} \to \mathbf{x}_j$, that is, the midpoint of the element. This limit gives

$$\eq{
\mathbb{G}(\mathbf{x}_j, \mathbf{x}_j) = \frac{a_j}{\pi} (\ln a_j - 1).
}$$

### Gauss-Legendre quadrature
The [Gauss-Legendre quadrature][GLquad] is a numerical integration method. For integrating a function of one variable $f(x)$ over the interval $[-1, 1]$, the rule takes the form:

$$\eq{
\int_{-1}^{1} f(x) \,dx \approx \sum_{i=1}^{n} w_i f(x_i),
}$$

where
* $n$, also called the quadrature order, is the number of sample points used
* $w_i$ are quadrature weights
* $x_i$ are the roots of the $n$-th Legendre polynomial

The abcissas $x_i$ and corresponding weights $w_i$ can be easily found in books tables or on internet pages, for example [here][xiwi]. For $n = 4$, which is the order used in this article, they are:

<center>

| $i$ |        $x_i$        |        $w_i$       |
|:---:|:-------------------:|:------------------:|
|  1  | -0.3399810435848563 | 0.6521451548625461 | 
|  2  |  0.3399810435848563 | 0.6521451548625461 | 
|  3  | -0.8611363115940526 | 0.3478548451374538 | 
|  4  |  0.8611363115940526 | 0.3478548451374538 | 

</center>

Applying the appropriate change of variables to change the integration interval to $[-1, 1]$, $(9)$ and $(10)$ can be rewritten as

$$\eq{
\mathbb{G}(\mathbf{x}_j, \boldsymbol{\xi}) = 
a_j \int_{-1}^{1} G(a_j \eta - \xi_x, \xi_y) \,d\eta \approx a_j \sum_{i=1}^{n} w_i G(a_j x_i - \xi_x, \xi_y),
}$$

$$\eq{
\mathbb{Q}(\mathbf{x}_j, \boldsymbol{\xi}) = 
a_j \int_{-1}^{1} Q(a_j \eta - \xi_x, \xi_y) \,d\eta \approx a_j \sum_{i=1}^{n} w_i Q(a_j x_i - \xi_x, \xi_y).
}$$

### Scipy.integrate.quad

## Results

{{< figure src="images/test_points.svg" alt="Test points" align="center" >}}

{{< figure src="images/GQ_abs_error.svg" alt="G and Q absolute error" align="center" >}}

{{< figure src="images/eval_time.svg" alt="Evaluation time" align="center" >}}

## Conclusion


## References
1. Steven L. Crouch, Sofia G. Mogilevskaya. 2024. A First Course in Boundary Element Methods. Springer, Cham, Switzerland.

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[2_bem_python]: ../2_bem_python
[spquad]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
[GLquad]: https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
[fcbem]: https://search.worldcat.org/title/1450559181
[xiwi]: https://pomax.github.io/bezierinfo/legendre-gauss.html
