---
date: '2025-09-28'
draft: false
title: 'Introducing Finite Element Method concepts with Julia'
author: 'Rodrigo Castro'
summary: 'The Julia programming Language is used to introduce basic concepts in Finite Element Methods by solving a two-dimensional truss problem.'
tags: ['Julia', 'Finite Element Method']
---

## Introduction
The objective of this article is to present some introductory concepts in Finite Element Methods by solving a two-dimensional truss problem. The inspiration of this post is the first chapter of the book *A First Course in Finite Elements* by *Jacob Fish and Ted Belytschko*. The [Julia] programming language is used to implement the algorithm suggested by the book, and the results are complemented with those obtained by [Ftool], a software for structural analysis.

## Methods
This topic starts with the construction of the element and truss stiffness matrices and its implementation in Julia. Next, the treatment of boundary conditions and the solution for nodal displacements and forces are presented. Finally, The truss to be evaluated as an example is defined.

### Bar element
The following image displays the configuration for a bar element and its two nodes, their displacements and the forces acting on them.

{{< figure src="images/bar_element.svg" alt="element" align="center" >}}

The local stiffness matrix $\mathbf{K}^e$ relates nodal forces and nodal displacements for the element $e$ in the global coordinate system $xy$ as

$$\eq{
\mathbf{F}^e = \mathbf{K}^e \mathbf{d}^e,
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

In the expressions above for the bar element, $E^e$ is the Young's Modulus, $A^e$ is the cross-section area and $l^e$ is the element length.

The concepts above are implemented in the following Julia composite type:

{{< dropdown_file title="Julia Element type" src="type_element.jl" fmt="julia" >}}

### Truss
A truss is a collection of interconnected bar elements and its stiffness matrix is constructed by combining the contributions of each element. In order to build the truss stiffness matrix $\mathbf{K}$, each element stiffness matrix $\mathbf{K}^e$ is first transformed from its local representation displayed in equation $(2)$ to a global representation, and the elements contributions are added together as

$$\eq{
\mathbf{K} = \sum_{e=1}^{n_{el}} \mathbf{L}^{e \mathrm{T}} \mathbf{K}^e \mathbf{L}^e.
}$$

$\mathbf{L}^e$ is called *gather* matrix, it has 4 rows (2 for each element node) and $2 n_\mathrm{nodes}$ columns, where $n_\mathrm{nodes}$ is the total number of nodes in the truss. This matrix is almost entirely composed of zeros, being ones only the entries that map the element nodes to their global numbering.

The following implementation in Julia shows the construction of the truss stiffness matrix in a composite type:

{{< dropdown_file title="Julia Truss type" src="type_truss.jl" fmt="julia" >}}

### Boundary value problem
Some of the truss nodes are fixed (*essential* nodes) and some are *free* to displace. In addition to that, external forces are applied to some of the nodes and these will be the loads that make the truss deform. To aid the solution process, the equation that relates forces and displacements is arranged in the following form:

$$\eq{
\begin{bmatrix}
\mathbf{K}^\mathrm{E} & \mathbf{K}^\mathrm{EF} \\
\mathbf{K}^\mathrm{FE} & \mathbf{K}^\mathrm{F} \\
\end{bmatrix}
\begin{bmatrix}
\mathbf{d}^\mathrm{E} \\
\mathbf{d}^\mathrm{F}
\end{bmatrix} = 
\begin{bmatrix}
\mathbf{0} \\
\mathbf{f}^\mathrm{F}
\end{bmatrix} + 
\begin{bmatrix}
\mathbf{r}^\mathrm{E} \\
\mathbf{0}
\end{bmatrix},
}$$

where the superscript $\mathrm{E}$ stands for *essential* and $\mathrm{F}$ stands for *free*. In this article, essential nodes are considered fixed, therefore $\mathbf{d}^\mathrm{E} = \mathbf{0}$. So, the displacement of free nodes $\mathbf{d}^\mathrm{F}$ can be computed as

$$\eq{
\mathbf{K}^\mathrm{F} \mathbf{d}^\mathrm{F} = \mathbf{f}^\mathrm{F}.
}$$

After solving for $\mathbf{d}^\mathrm{F}$, the reactions at the essential nodes can be computed from

$$\eq{
\mathbf{r}^\mathrm{E} = \mathbf{K}^\mathrm{EF} \mathbf{d}^\mathrm{F},
}$$

and the tensile stress can be computed from $(3)$.

All these steps are implemented in Julia through the following function:

{{< dropdown_file title="Julia solve function" src="function_solve.jl" fmt="julia" >}}

### Application problem
The code implemented in Julia will be tested against the truss presented in following figure. I consists of 11 steel bar elements, all with circular cross-section with 2 cm of diameter. The complete truss has 3 m in length and 1 m in height. The elements are connected by 7 nodes, where nodes 1 and 4 are fixed, but free to rotate. Vertical loads of -0.5 kN and -1.0 kN are concentrated at nodes 2 and 3, respectively.

{{< figure src="images/ftool_truss.svg" alt="truss" align="center" >}}

In Julia, the truss is defined by the following code:

{{< dropdown_file title="Julia truss definition" src="application.jl" fmt="julia" >}}

In the code above, the nodes are defined by their *x* and *y* coordinates in the *nodes* matrix. Elements are defined by each line of the *connectivity* matrix, where the two columns represent the first and second nodes, respectively.

## Results
The truss problem defined in the last topic is solved and the solution is complemented with results obtained by [Ftool], a software for structural analysis. Their results will not be compared quantitatively side by side because they are actually identical. So, the following images, displaying results from the code implemented in Julia and from Ftool, represent the solution obtained throught the Finite Element Method.

The first two images displays the truss in deformed condition, scaled by a factor of 1000. The first one computed with Julia and plotted with [Makie], the second evaluated by Ftool.

{{< figure src="images/truss_deformation.svg" alt="truss deformation" align="center" >}}

{{< figure src="images/ftool_deformation.svg" alt="ftool truss deformation" align="center" >}}

The next two images shows the actual displacements in *x* and *y*, respectively.

{{< figure src="images/ftool_dx.svg" alt="ftool dx" align="center" >}}

{{< figure src="images/ftool_dy.svg" alt="ftool dx" align="center" >}}

Ftool does not compute the tensile stress directly, but rather the elements axial forces. These are displayed by the following two images, the first computed with Julia and plotted with [Makie], the second evaluated by Ftool.

{{< figure src="images/axial_load.svg" alt="element load" align="center" >}}

{{< figure src="images/ftool_load.svg" alt="ftool element load" align="center" >}}

## Conclusion
This post presents an interesting and fun way of introducing Finite Element concepts. The truss problem is easily implemented in any programming language and showcases why the Finite Element Methods is so dominant in structural analysis. However, this introduction does not tell much what is FEM and what it can be used for. Future articles will hopefully make this clearer.

## References
1. Jacob Fish and Ted Belytschko. 2007. A First Course in Finite Elements. Wiley, Chichester, UK.

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[Julia]: https://julialang.org/
[Ftool]: https://www.ftool.com.br/Ftool/
[Makie]: https://docs.makie.org/
