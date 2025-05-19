---
date: '2025-05-19'
draft: false
title: '2D constant boundary element influence coefficients'
author: 'Rodrigo Castro'
summary: 'Computation of influence coefficients for 2D constant boundary elements.'
tags: ['Boundary Element Method']
---

## Introduction
This post presents a first improvement on the [previous article][2_bem_python], specially on the calculation of the influence coefficients for the two dimensional boundary element method.

The ojective is to study how the [`scipy.integrate.quad`][spquad] function and the [Gauss-Legendre quadrature][GLquad] compare to the analytical evaluation of the integrals that represent the influence coefficients.

## Methods
The [previous post][2_bem_python] introduced the boundary integral equation for the Laplace's equation. The equation is presented below with a similar nomenclature to that used by *Crouch and Mogilevskaya*[<sup>[1]</sup>](#references) in their [book][fcbem], which will be the main reference for this work.

$$\eq{
\frac{1}{2} \phi(\boldsymbol{\xi}) = 
\int_{S} \phi(\mathbf{x}) Q(\mathbf{x}, \boldsymbol{\xi}) \,dS(\mathbf{x}) -
\int_{S} q(\mathbf{x}) G(\mathbf{x}, \boldsymbol{\xi}) \,dS(\mathbf{x}),
}$$

in which,

$$\eq{
q(\mathbf{x}) = \frac{\partial \phi(\mathbf{x})}{\partial \mathbf{n}(\mathbf{x})},
}$$

$$\eq{
Q(\mathbf{x}, \boldsymbol{\xi}) = \frac{\partial G(\mathbf{x}, \boldsymbol{\xi})}{\partial \mathbf{n}(\mathbf{x})}.
}$$

This formula expresses the potential $\phi(\boldsymbol{\xi})$ at a point $\boldsymbol{\xi}$ within an homogeneous region $V$ in terms of integrals of the potential $\phi(\mathbf{x})$ and its normal derivative $q(\mathbf{x})$ over the boundary $S$ of this region.

$G(\mathbf{x}, \boldsymbol{\xi})$ is the fundamental solution to Laplace's equation that gives the potential at the sample point $\mathbf{x}$ due to a source point at $\boldsymbol{\xi}$. The function $Q(\mathbf{x}, \boldsymbol{\xi})$ is the flux across the boundary $S$ at point $\mathbf{x}$ associated with the point $\boldsymbol{\xi}$. In two dimensions, the fundamental solution $G$ is: 

$$\eq{
G(\mathbf{x}, \boldsymbol{\xi}) = \frac{1}{2 \pi} \ln |\mathbf{x} - \boldsymbol{\xi}|.
}$$

The use of the <abbr title="Boundary Element Method">BEM</abbr> code that was implemented in the [last article][2_bem_python] requires the discretization of the boundary $S$ into $N$ straight elements with constants $\phi$ and $q$ over each element. Then, the boundary integral equation in $(1)$ can be approximated by

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

The following subtopics present how to evaluate the integrals $(6)$ and $(7)$ by three methods. The first method is analytical and the other two are numerical, using the Gauss-Legendre quadrature and the `scipy.integrate.quad` function.

### Analytical integration
Integrals $(6)$ and $(7)$ are evaluated over a straight boundary element of length $2a_j$, shown as a red color line in the following figure. The element local system has the $x$ axis colinear with the element and the $y$ axis normal to it. The blue dot represents the source point $\boldsymbol{\xi}$, which has coordinates $(\xi_x, \xi_y)$ in the element local system. The other displayed paramaters will simplify the resulting analytical expression and are given in equations $(8)$ and $(9)$.

{{< figure src="images/intGQ.svg" alt="analytical G Q parameters" align="center">}}

$$\eq{
r_1 = \sqrt{\left(\xi_x - a_j \right)^2 + \xi_y^2}, \quad r_2 = \sqrt{\left(\xi_x + a_j \right)^2 + \xi_y^2},
}$$

$$\eq{
\theta_1 = \arctan \frac{\xi_y}{\xi_x - a_j}, \quad \theta_2 = \arctan \frac{\xi_y}{\xi_x + a_j}.
}$$

$\theta_1$ and $\theta_2$ can assume values between $-\pi$ and $\pi$.

That being defined, integrals $(6)$ and $(7)$ can be rewritten as

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
G(x - \xi_x, \xi_y) = 
\frac{1}{2 \pi} \ln \sqrt{\left( x - \xi_x \right)^2 + \xi_y^2},
}$$

$$\eq{
Q(x - \xi_x, \xi_y) = 
-\frac{1}{2 \pi} \frac{\xi_y}{\left( x - \xi_x \right)^2 + \xi_y^2}.
}$$

Expressions $(10)$ and $(11)$ can be analytically evaluated to

$$\eq{
\mathbb{G}(\mathbf{x}_j, \boldsymbol{\xi}) = 
\frac{1}{2 \pi} \left[ \xi_y (\theta_1 - \theta_2) - (\xi_x - a_j) \ln r_1 + (\xi_x + a_j) \ln r_2 - 2a_j \right],
}$$

$$\eq{
\mathbb{Q}(\mathbf{x}_j, \boldsymbol{\xi}) = 
-\frac{1}{2 \pi} (\theta_1 - \theta_2).
}$$

For any $\boldsymbol{\xi}$ colinear with the element, $\mathbb{Q}$ is evaluated to

$$\eq{
\mathbb{Q}(\mathbf{x}_j, \boldsymbol{\xi}) = 0.
}$$

The value of $\mathbb{G}$ for $\boldsymbol{\xi}$ coincident with the element's endpoints $(\pm a_j, 0)$ is obtained considering the limit case of the expression $(14)$ as $\xi_y = 0$ and $\xi_x \to \pm a_j$. 

$$\eq{
\mathbb{G}(\mathbf{x}_j, (\pm a_j, 0)) = \frac{a_j}{\pi} [\ln(2 a_j) - 1].
}$$

Similarly, the diagonal (self-effect) terms in $(5)$ can be obtained from the limit case of $(14)$ as $\boldsymbol{\xi} \to \mathbf{x}_j$, that is, the midpoint of the element. This limit is

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

The abcissas $x_i$ and corresponding weights $w_i$ can be easily found in books' tables or on internet pages, for example [here][xiwi]. For $n = 4$, which is the order used in this work, they are:

<center>

| $i$ |        $x_i$        |        $w_i$       |
|:---:|:-------------------:|:------------------:|
|  1  | -0.3399810435848563 | 0.6521451548625461 | 
|  2  |  0.3399810435848563 | 0.6521451548625461 | 
|  3  | -0.8611363115940526 | 0.3478548451374538 | 
|  4  |  0.8611363115940526 | 0.3478548451374538 | 

</center>

Applying the appropriate change of variables to change the integration interval to $[-1, 1]$, $(10)$ and $(11)$ can be approximated by the following expressions:

$$\eq{
\mathbb{G}(\mathbf{x}_j, \boldsymbol{\xi}) = 
a_j \int_{-1}^{1} G(a_j \eta - \xi_x, \xi_y) \,d\eta \approx a_j \sum_{i=1}^{n} w_i G(a_j x_i - \xi_x, \xi_y),
}$$

$$\eq{
\mathbb{Q}(\mathbf{x}_j, \boldsymbol{\xi}) = 
a_j \int_{-1}^{1} Q(a_j \eta - \xi_x, \xi_y) \,d\eta \approx a_j \sum_{i=1}^{n} w_i Q(a_j x_i - \xi_x, \xi_y).
}$$

### Scipy.integrate.quad
[`quad`][spquad] is a general purpose integration function that is part of [SciPy], a Python package that provides algorithms for scientific computing. An example of its use is given in the following snippet, which computes $\mathbb{G}$ and $\mathbb{Q}$ for specific values of $a_j$ and $\boldsymbol{\xi}$.

```python
import numpy as np
from scipy.integrate import quad

aj = 0.5
cx, cy = [1.0, 1.0]  # ξ

intG = lambda x: np.log((cx - x)**2 + cy**2)
intQ = lambda x: cy / ((cx - x)**2 + cy**2)

singular_point = [cx] if cy == 0 else None

G = 0.25 / np.pi * quad(intG, -aj, aj, points=singular_point)[0]
Q = -0.5 / np.pi * quad(intQ, -aj, aj, points=singular_point)[0]
```

`quad` accepts a sequence of breakpoints (`points`) in the integration interval where local difficulties of the integrand may occur. It evaluates the integral iteratively until the absolute or relative error gets below a given value, which by default are set to `epsabs=1.49e-8` and `epsrel=1.49e-8`. Also, it returns an array with two elements, the integral and an estimate of the absolute error.

## Results
To compare the three integration methods, a Python code calculates the values of $\mathbb{G}$ and $\mathbb{Q}$ for 1000 sources $\boldsymbol{\xi}$ spread randomly around the element $\mathbf{x}$, as depicted in the figure below (the image does not show all the source points). The element has unitary length and the distances between the sources and the element are bounded by $|\mathbf{x} - \boldsymbol{\xi}| ≤ 10$.

{{< figure src="images/test_points.svg" alt="Test points" align="center" >}}

The next image presents the absolute error of `quad` and Gauss-Legendre quadrature of order 4 compared to the analytical results. It can be seen from the two bottom charts that the Gauss-Legendre quadrature loses precision the closer the source gets to the element. The same can be seen from `quad`, but less intensily due to the adaptive behavior of the function.

{{< figure src="images/GQ_abs_error.svg" alt="G and Q absolute error" align="center" >}}

The following picture presents the amount of seconds it take to evalute $\mathbb{G}$ and $\mathbb{Q}$ as a function of the distance between the source and the element. Again, the adaptive behavior of `quad` explains the longer time when the source is close to the element. On average, `quad` takes 3 times longer than the analytical method, while the Gauss-Legendre quadrature of order 4 takes 1.5 times longer than the analytical approach.

{{< figure src="images/eval_time.svg" alt="Evaluation time" align="center" >}}

## Conclusion
This study encourages the use of the analytical approach to improve the Boundary Element Method code implemented in the [previous post][2_bem_python], as it is the faster and the most precise. The Gauss-Legendre quadrature will be kept as an alternative and the probable main method for future tasks that require complex Green's functions or more sofisticated boundary elements, where the analytic integration method is time consuming or inexistent. However, the Gauss-Legendre quadrature implementation still requires polishment, inspired by `quad`:

* Better approximations for the influence coefficients when $\boldsymbol{\xi}$ is close to $\mathbf{x}$;
* Handling of $G$ singularities in the integration domain;
* Make computations faster by implementing it in a high performance programming language, like Fortran ([See how Fortran and Python can be integrated][1_f2py]).

## References
1. Steven L. Crouch, Sofia G. Mogilevskaya. 2024. A First Course in Boundary Element Methods. Springer, Cham, Switzerland.

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[1_f2py]: ../1_f2py_fortran_python
[2_bem_python]: ../2_bem_python
[spquad]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
[GLquad]: https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
[fcbem]: https://search.worldcat.org/title/1450559181
[xiwi]: https://pomax.github.io/bezierinfo/legendre-gauss.html
[SciPy]: https://scipy.org/
