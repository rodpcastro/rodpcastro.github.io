---
date: '2025-09-27'
draft: true
title: 'Introducing Finite Element Method concepts with Julia'
author: 'Rodrigo Castro'
summary: 'The Julia programming Language is used to introduce basic concepts in Finite Element Methods by solving a two-dimensional truss problem.'
tags: ['Julia', 'Finite Element Method', 'Ftool']
---

## Introduction
The objective of this article is to present some introductory concepts in Finite Element Methods by solving a two-dimensional truss problem. The inspiration of this post is the first chapter of the book *A First Course in Finite Elements* by *Jacob Fish and Ted Belytschko*. The [Julia] programming language is used to implement the algorithm suggested by the book, and the results are compared with those obtained by [Ftool], a software for structural analysis of trusses.

## Methods
This topic starts with the construction of the truss stiffness matrix and its implementation in Julia. Next, the treatment of boundary conditions and the solutions for nodal displacements and forces are presented. Finally, The truss to be evaluated is defined.

### Bar element
The following image displays a bar element configuration with two nodes, their displacements and the forces acting on them.

{{< figure src="images/bar_element.svg" alt="element" align="center" >}}

The local stiffness matrix $\mathbf{K}^e$ relates nodal forces and nodal displacements for the element $e$ in the global coordinate system $xy$ as

$$\eq{
\mathbf{F}^e = \mathbf{K}^e \mathbf{d}^e
}$$

where

$$\eq{
\mathbf{K}^e = \frac{E^e A^e}{l^e}
\begin{bmatrix}
c^2 & s c & -c^2 & -s c \\
s c & s^2 & -s c & -s^2 \\
-c^2 & -s c & c^2 & s c \\
-s c & -s^2 & s c & s^2 \\
\end{bmatrix},
}$$

and

$$
s^2 = \sin^2\phi^e, \quad c^2 = \cos^2\phi^e, \quad s c = \sin\phi^e \cos\phi^e.
$$

The tensile stress is computed, after solving for $\mathbf{d}^e$, as

$$\eq{
\sigma^e = \frac{E^e}{l^e} 
\begin{bmatrix}
-c & -s & c & s
\end{bmatrix}
\mathbf{d}^e.
}$$

In the expressions above for the bar element, $E^e$ is the Young's Modulus, $A^e$ is the cross-section area and $l^e$ is the length.

The concepts above are implemented in the following Julia composite type:

{{< dropdown_file title="Julia Element type" src="type_element.jl" fmt="julia" >}}

### Truss

{{< dropdown_file title="Julia Truss type" src="type_truss.jl" fmt="julia" >}}

### Boundary value problem

{{< dropdown_file title="Julia solve function" src="function_solve.jl" fmt="julia" >}}

### Application problem

{{< figure src="images/ftool_truss.svg" alt="truss" align="center" >}}

## Results


## Conclusion


## References
1. Jacob Fish and Ted Belytschko. 2007. A First Course in Finite Elements. Wiley, Chichester, UK.


## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[Julia]: https://julialang.org/
[Ftool]: https://www.ftool.com.br/Ftool/
