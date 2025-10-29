---
date: '2025-10-29'
draft: false
title: 'Finding natural frequencies and modes of a vibrating beam with the Rayleigh-Ritz method'
author: 'Rodrigo Castro'
summary: 'The Rayleigh-Ritz method is used to find the natural frequencies and modes of a propped cantilever beam.'
tags: ['Variational Calculus', Rayleigh-Ritz', 'Structural Dynamics']
---

## Introduction
This post is the first on the topic of Variational Calculus. As an application of it, I decided to use the Rayleigh-Ritz method to find the natural frequencies and modes of a propped cantilever beam.

## Methods
The Rayleigh-Ritz method is explained first, followed by it's application in the case a propped cantilever beam.

### Rayleigh-Ritz method
The obtainment of the natural frequencies and modes of an oscillating system by the Rayleigh-Ritz method is fundamented on the Rayleigh's principle, which is embedded in the following quote:

> Equating the maximum total potential energy associated with vibration to the maximum kinetic energy associated with vibration results in an upperbound estimate of the fundamental natural frequency, provided the dynamic displacement forms assumed are admissible. A displacement function is admissible if it does not violate any geometric constraints and can represent the displaced form of the system without any discontinuity.
> - <cite>Sinniah Ilanko</cite>

Although the above paragraph only speaks about the fundamental natural frequency, which is the lowest one, the Rayleigh-Ritz method employs the same principle to derive upper bounds to any natural frequency of an oscillating system. The method is described as follows.

Let $f$ be the maximum displacement for a one dimensional continuous system, expressed as a series of products of undetermined weighting coefficients and admissible displacement functions $\phi_i$.

$$\eq{
f = \sum G_i \phi_i,
}$$

The maximum potential energy $V_m$ and maximum kinetic energy $T_m$ can also be expressed as functions of $G_i$ and $\phi_i$. If we write

$$\eq{
T_m = \omega^2 \psi_m,
}$$

where $\psi_m$ is a kinetic energy function, then according to the Rayleigh's principle,

$$\eq{
\omega^2 = \frac{V_m}{\psi_m}.
}$$

To obtain the best possible estimate for other natural frequencies other than the fundamental, we derive equation $(3)$ with respect to $G_i$ and make it equal to zero, resulting in the equation

$$\eq{
\frac{\partial V_m}{\partial G_i} - \omega^2 \frac{\partial \psi_m}{\partial G_i} = 0.
}$$

If the number of terms in the series for $f$ is $n$, there will be $n$ such equations, which for linear systems can be written in matrix form as follows

$$\eq{
[K]\{G\} - \omega^2 [M] \{G\} = \{0\}.
}$$

The coefficients of $[K]$ and $[M]$ are given by

$$\eq{
K_{ij} = \frac{\partial^2 V_m}{\partial G_i G_j},
}$$

$$\eq{
M_{ij} = \frac{\partial^2 \psi_m}{\partial G_i G_j}.
}$$

The solution to the eigenvalue problem $(5)$ gives $n$ natural frequencies (eigenvalues) and $n$ corresponding modes (eigenvectors).

### Expressions for a propped cantilever beam
The Rayleigh-Ritz method will be applied for a propped cantilever beam, displayed in the following figure, where $L$ is beam's length, $E$ is the Young's modulus, $I$ is the area moment of inertia and $m$ is the linear mass.

{{< figure src="images/beam.svg" alt="Beam" align="center" >}}

The geometric constraints of this problem are

$$\eq{
\begin{split}
f(0) = 0, \\
f'(0) = 0, \\
f(L) = 0,
\end{split}
}$$

and the admissible displacement functions $\phi_i$ that satisfy these constraints are given by

$$\eq{
\phi_i = \bar{x}^{i+2} (1-\bar{x}), \quad i = 0, 1, \ldots, n
}$$

where $\bar{x} = x / L$. $V_m$ and $\psi_m$ for the beam are expressed as

$$\eq{
V_m = \frac{E I}{2 L^3} \int_{0}^{1} (\bar{f}''(\bar{x}))^2 \, d\bar{x},
}$$

$$\eq{
\psi_m = \frac{m L}{2} \int_{0}^{1} (\bar{f}(\bar{x}))^2 \, d\bar{x}.
}$$

Substituting $(10)$ and $(11)$ into equation $(5)$, and after some manipulations we get

$$\eq{
[\bar{K}]\{G\} - \lambda^4 [\bar{M}] \{G\} = \{0\},
}$$

where $\lambda$ is a dimensionless frequency parameter defined by

$$\eq{
\lambda = \left(\frac{\omega^2 m L^4}{E I}\right)^{1/4},
}$$

and the elements of the matrices $[\bar{K}]$ and $[\bar{M}]$ are

$$\eq{
\bar{K}_{ij} = 2 \int_{0}^{1} \phi_i'' \phi_j'' \, d\bar{x},
}$$

$$\eq{
\bar{M}_{ij} = 2 \int_{0}^{1} \phi_i \phi_j \, d\bar{x},
}$$

which were evaluated algebraically with [`sympy`][sympy] to give

$$\eq{
\bar{K}_{ij} = \frac{4 \left(i + 2\right) \left(j + 2\right) \left(i^{2} + 3 i j + 4 i + j^{2} + 4 j + 3\right)}{\left(i + j + 1\right) \left(i + j + 2\right) \left(i + j + 3\right)},
}$$

$$\eq{
\bar{M}_{ij} = \frac{4}{\left(i + j + 5\right) \left(i + j + 6\right) \left(i + j + 7\right)}.
}$$

## Results
Natural frequencies and modes obtained with the Rayleigh-Ritz method are compared to the analytical solution, extracted from chapter 4 of the book by *Blevins* (2016). The first figure below shows the first four natural frequency parameters as a function of the number of terms used in equation $(1)$. The second figure displays the first four natural modes for eight terms used.

{{< figure src="images/freqparams.svg" alt="frequency parameters" align="center" >}}

{{< figure src="images/modes.svg" alt="natural modes" align="center" >}}

As demonstrated by the two figures above, the Rayleigh-Ritz method with eight terms provide great approximation for the propped cantilever beam natural frequencies and modes.

## Conclusion
The Rayleigh-Ritz method is a simple but useful application on the topic of Variational Calculus. Future posts will hopefully bring more examples on which this topic is important to understand problems in physics and engineering.

## References
1. Sinniah Ilanko. 2014. The Rayleigh-Ritz Method for Structural Analysis. Wiley-ISTE.
2. Robert D. Blevins. 2016. Natural Frequency of Beams. In Formulas for Dynamics, Acoustics and Vibration. John Wiley & Sons, 134–202.

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[sympy]: https://www.sympy.org/
