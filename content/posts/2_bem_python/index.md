---
date: '2025-05-05'
draft: false
title: 'Introducing the Boundary Element Method with Python'
author: 'Rodrigo Castro'
summary: 'Python implementation of the 2D Boundary Element Method.'
tags: ['Boundary Element Method', 'Python']
---

## Introduction
The [Boundary Element Method] (BEM) is a numerical computational technique for solving partial differential equations in physics and engineering. Its main advantage over other numerical methods, such as the Finite Element Method, is its ability to reduce the problem to the boundary, providing a solution with n-1 dimensions to a problem of n dimensions. This makes it highly efficient and often the preferred method for many applications in fluid mechanics, acoustics and electromagnetics.

This post represents my introduction to the 2D Boundary Element Method by reproducing the [article] of *Keng-Cheng Ang*. Unlike the original, which uses MATLAB, I implement the method in Python.

## Methods

### Laplace's equation
The  aim is to solve the [Laplace's equation] in the 2D region \(R\), subject to Dirichlet and Neumann boundary conditions in \(C_\alpha\) and \(C_\beta\), respectively:

$$ \nabla^2 u = \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = 0, \quad \text{for } (x, y) \in R, \tag{1} $$

$$ u = f(x,y) \quad \text{for } (x, y) \in C_\alpha, \tag{2} $$

$$ q = \frac{\partial u}{\partial n} = g(x, y) \quad \text{for } (x, y) \in C_\beta. \tag{3} $$

The fundamental solution of Laplace's equation in two dimensions is:

$$ u^* = \frac{1}{4\pi} \ln \left( (x - \xi)^2 + (y - \eta)^2 \right). \tag{4} $$

Through application of [Green's Second Identity][greensids], it's possible to write the boundary integral equation for the Laplace's equation as:

$$ \lambda u(\xi, \eta) = \int_{C} \left( u \frac{\partial u^*}{\partial n} - u^* \frac{\partial u}{\partial n} \right) \,ds, \tag{5} $$

where \(C = C_\alpha \cup C_\beta\) and

$$
\lambda = \begin{cases}
  0 & \text{if } (\xi, \eta) \notin R \cup C,\\
  \frac{1}{2} & \text{if } (\xi, \eta) \text{ lies on a smooth part of } C,\\
  1 & \text{if } (\xi, \eta) \in R.
\end{cases}
\tag{6}
$$

With BEM, it is possible to find approximations to all the unkowns on the right-hand side of \((5)\). This is described next.

### The Boundary Element Method
The first step is approximating the boundary \(C\) by a polygon with \(N\) straight line segments:

$$ C \approx C_1 \cup C_2 \cup C_3 \cup \ldots \cup C_N \tag{7} $$

A closed boundary with \(N\) segments is defined by the \((x, y)\) coordinates of \(N+1\) nodes, with the last node being coincident with the first, that is \((x_{N+1}, y_{N+1}) = (x_1, y_1)\). A segment \(C_k\) is defined by a straight line joining \((x_k, y_k)\) to \((x_{k+1}, y_{k+1})\). These segments are known as boundary elements.

Our next simplification is taking the values of \(u\) and \(\partial u / \partial n\) in each element to be constant and equal to their values at the element midpoint \((\bar{x}_k, \bar{y}_k)\). That is,

$$ u \approx \bar{u}_k, \quad \text{and} \quad \frac{\partial u}{\partial n} = \bar{q}_k \quad \text{for } (x, y) \in C_k. \tag{8} $$

Hence, equation \((5)\) becomes:

$$ \lambda u(\xi, \eta) = \sum_{k=1}^{N} \left( \bar{u}_k G_k(\xi, \eta) - \bar{q}_k F_k(\xi, \eta) \right) \tag{9} $$

where

$$ F_k(\xi, \eta) = \int_{C_k} u^*(\xi, \eta) \,ds, \tag{10} $$

$$ G_k(\xi, \eta) = \int_{C_k} \frac{\partial u^*(\xi, \eta)}{\partial n} \,ds. \tag{11} $$

In a properly posed boundary value problem, either \(\bar{u}_k\) or \(\bar{q}_k\) (but not both) is known on any given point of the boundary. Equation \((9)\) can be used to find approximations to the unkown boundary values. To do so, let \((\xi, \eta)\) be the midpoint of element \(C_k\), which gives \(N\) equations with \(N\) unkowns:

$$
\frac{1}{2} \bar{u}_m = \sum_{k=1}^{N}
\left( \bar{u}_k G_k(\bar{x}_m, \bar{y}_m) - \bar{q}_k F_k(\bar{x}_m, \bar{y}_m) \right) \quad \text{for } m = 1, 2, \ldots, N,
\tag{12}
$$

where

$$ u^*_m = \frac{1}{4\pi} \ln \left( (x - \bar{x}_m)^2 + (y - \bar{y}_m)^2 \right). \tag{13} $$

Equation \((12)\) can be conveniently written in the matrix form \(A z = b\) by placing unkowns \(\bar{u}\) and \(\bar{q}\) in \(z\). As a consequence, the elements of \(A\) and \(b\) are given by

$$
A_{mk} = \begin{cases}
  -F_k(\bar{x}_m, \bar{y}_m) & \text{if } u \text{ given over } C_k,\\
  G_k(\bar{x}_m, \bar{y}_m) & \text{if } q \text{ given over } C_k \text{ and } k \ne m,\\
  G_k(\bar{x}_m, \bar{y}_m) - \frac{1}{2} & \text{if } q \text{ given over } C_k \text{ and } k = m.
\end{cases}
\tag{14}
$$

$$
b_{mk} = \begin{cases}
  \bar{q}_k F_k(\bar{x}_m, \bar{y}_m) & \text{if } q \text{ given over } C_k,\\
  \bar{u}_k G_k(\bar{x}_m, \bar{y}_m) & \text{if } u \text{ given over } C_k \text{ and } k \ne m,\\
  \bar{u}_k \left( -G_k(\bar{x}_m, \bar{y}_m) + \frac{1}{2} \right) & \text{if } u \text{ given over } C_k \text{ and } k = m.
\end{cases}
\tag{15}
$$

Once \(A\) and \(b\) are formed, the system can be solved for \(z\), bearing in mind that \(z_k = \bar{u}_k\) if \(q\) is given over \(C_k\), and \(z_k = \bar{q}_k\) if \(u\) is given over \(C_k\).

With all values of \(\bar{u}_k\) and \(\bar{q}_k\), Equation \((9)\) can be used to find the value of \(u\) at any interior point in the domain \(R\).

### Python implementation

The Python function below defines the boundary shape and the boundary conditions' types and values for the application used in the original article and also in this post. This application is better described in the section [Example](#example). Here, we simply explain what are each one the variables returned by this function represent.

{{< dropdown_file title="Boundary definition" src="files/boundary_definition.py" fmt="python" >}}

In the function above,

* **xb** and **yb** are vector arrays containing \(x\) and \(y\) coordinates of \(N+1\) elements' nodes.
* **bt** is a vector array containing the boundary coditions' types of \(N\) elements. \(0\) and \(1\) represent the Dirichlet and Neumann boundary conditions, respectively.
* **bv** is a vector array containing the boundary conditions' values of \(N\) elements.



### Example
?

## Results

{{< figure src="analytical_bem_comparison.svg" alt="analytical bem comparison" align="center" >}}

{{< figure src="abserr.svg" alt="bem absolute error" align="center" >}}

## Conclusion


## References
1. Keng-Cheng Ang. 2008. Introducing the boundary element method with MATLAB. International Journal of Mathematical Education in Science and Technology 39, 4 (Jun. 2008), 505–19. https://doi.org/10.1080/00207390701722676

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[Boundary Element Method]: https://en.wikipedia.org/wiki/Boundary_element_method 
[article]: https://doi.org/10.1080/00207390701722676
[Laplace's equation]: https://en.wikipedia.org/wiki/Laplace%27s_equation
[greensids]: https://en.wikipedia.org/wiki/Green%27s_identities
