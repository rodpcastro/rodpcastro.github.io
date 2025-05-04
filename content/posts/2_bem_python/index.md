---
date: '2025-05-04'
draft: false
title: 'Introducing the Boundary Element Method with Python'
author: 'Rodrigo Castro'
summary: 'Python implementation of the 2D Boundary Element Method.'
tags: ['Boundary Element Method', 'Python']
---

## Introduction
The Boundary Element Method (BEM) is a numerical computational technique for solving partial differential equations in physics and engineering. Its main advantage over other numerical methods, such as the Finite Element Method, is its ability to reduce the problem to the boundary, providing a solution with n-1 dimensions to a problem of n dimensions. This makes it highly efficient and often the preferred method for many applications in fluid mechanics, acoustics and electromagnetics.

This post represents my introduction to the 2D Boundary Element Method by reproducing the [article] of *Keng-Cheng Ang*. Unlike the original, which uses MATLAB, I implement the method in Python.

## Methods


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
[article]: https://doi.org/10.1080/00207390701722676
