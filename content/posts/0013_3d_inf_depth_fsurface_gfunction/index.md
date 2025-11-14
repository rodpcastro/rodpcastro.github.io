---
date: '2025-10-13'
draft: false
title: 'Algorithms for the 3D infinite-depth free-surface Green function'
author: 'Rodrigo Castro'
summary: 'Series expansions for the three-dimensional infinite-depth free-surface Green function and its derivatives.'
tags: ['Green Function', 'Potential Flow', 'Python']
---

## Introduction
This post is the reproduction of an article by *Shan & Wu (2018)*, in which they derive series expansions for the three-dimensional infinite-depth free-surface Green function. The goal of this post is to set the initial steps for the future project of solving the floating body problem in three dimensions.

## Methods
The infinite-depth free-surface Green function represents the spatial component of a velocity potential induced at a field point $p(x, y, z)$ by a pulsating source point $q(\xi, \eta, \zeta)$. Being $\bar{q}(\xi, \eta, -\zeta)$ the image of $q$ relative to the free-surface, we define the following quantities:

$$
\begin{split}
r = \sqrt{(x-\xi)^2 + (y-\eta)^2} &, \quad Z = z+\zeta, \\
R_{pq} = \sqrt{r^2 + (z-\zeta)^2} &, \quad R_{p\bar{q}} = \sqrt{r^2 + Z^2}, \\
X = k_0 r, \quad Y = -k_0 Z &, \quad R = \sqrt{X^2 + Y^2}. \\
\end{split}
$$

where $k_0$ is the infinite-depth wave number. Between $p$ and $q$, $r$ is the horizontal distance and $R_{pq}$ is the Euclidian distance. Between $p$ and $\bar{q}$, $|Z|$ is the vertical distance and $R_{p\bar{q}}$ is the Euclidian distance. $X$ and $Y$ are non-negative dimensionless cylindrical coordinates, and $R$ is the dimensionless form of $R_{p\bar{q}}$.

The infinite-depth free-surface Green function $G_\infty(p,q)$ can be expressed as

$$\eq{
G_\infty(p,q) = \frac{1}{R_{pq}} + \frac{1}{R_{p\bar{q}}} + k_0 F(X,Y) \pm 2\mathrm{i}\pi k_0 e^{-Y} J_0(X),
}$$

where $J_0$ is the zero-order Bessel function of the first kind. The $(-)$ sign in the last expression is associated with the time component $e^{\mathrm{i} \omega t}$, while the $(+)$ sign corresponds to $e^{-\mathrm{i} \omega t}$, and $\omega$ is the pulsating source frequency.

The laborious evaluation of $G_\infty$ and its derivatives is translated to computing $F$ and its derivatives:

$$\eq{
\begin{split}
F = & -2 e^{-Y} \int_{0}^{Y} e^t (X^2+t^2)^{-\frac{1}{2}} \,dt \\
& -\pi e^{-Y} [H_0(X) + Y_0(X)],
\end{split}
}$$

$$\eq{
\begin{split}
\frac{\partial F}{\partial X} = & \phantom{-} 2 X e^{-Y} \int_{0}^{Y} e^t (X^2+t^2)^{-\frac{3}{2}} \,dt \\
& -2 e^{-Y} + \pi e^{-Y} [H_1(X) + Y_1(X)],
\end{split}
}$$

$$\eq{
\begin{split}
\frac{\partial^2 F}{\partial X^2} = & \phantom{-} 2 e^{-Y} \int_{0}^{Y} e^t (X^2+t^2)^{-\frac{5}{2}} (t^2-2X^2) \,dt \\
& +\frac{1}{3} e^{-Y} X + \frac{\pi}{2} e^{-Y} [H_0(X) + Y_0(X) - H_2(X) - Y_2(X)].
\end{split}
}$$

where $H_n$ is the $n$-th order Struve function and $Y_n$ is the $n$-th order Bessel function of the second kind. The derivatives with the respect to $Y$ are related to expressions above, so the focus is on the evaluation of $F$, $F_X$ and $F_{XX}$.

*Shan & Wu (2018)* derived series expansions for $F$ and its derivatives for four different regions of the first quadrant of the $XY$ plane. The following image depicts these regions $D_i$, for $i=1,\ldots,4$ and the number of terms $N_i$ that are used to approximate $F$ and its derivatives in each region.

{{< figure src="images/subdomains.svg" alt="Subdomains" align="center" >}}

The series expansions are given in the following subtopics.

### Series Expansion for $D_1$

$$\eq{
\begin{split}
F = &-2 e^{-Y} \mathrm{Ei}(Y) \\
& + 2 \sum_{n=1}^{\infty} \frac{1}{n!^2} \left(-\frac{X^2}{4}\right)^n
\left[\sum_{m=1}^{2n} \frac{(m-1)!}{Y^m} - e^{-Y} \mathrm{Ei}(Y)\right],
\end{split}
}$$

$$\eq{
\begin{split}
\frac{\partial F}{\partial X} = \frac{4}{X} \sum_{n=1}^{\infty} \frac{n}{n!^2} \left(-\frac{X^2}{4}\right)^n
\left[\sum_{m=1}^{2n} \frac{(m-1)!}{Y^m} - e^{-Y} \mathrm{Ei}(Y)\right],
\end{split}
}$$

$$\eq{
\begin{split}
\frac{\partial^2 F}{\partial X^2} = \frac{4}{X^2} \sum_{n=1}^{\infty} \frac{n(2n-1)}{n!^2} \left(-\frac{X^2}{4}\right)^n
\left[\sum_{m=1}^{2n} \frac{(m-1)!}{Y^m} - e^{-Y} \mathrm{Ei}(Y)\right].
\end{split}
}$$

$\mathrm{Ei}$ is the exponential integral.

### Series Expansion for $D_2$

$$\eq{
\begin{split}
F = &-\pi e^{-Y} Y_0(X) - \frac{2R}{X^2} \\
& + \frac{2R}{X^2} e^{-Y} \sum_{n=0}^{\infty} \frac{(n+1)X^n}{n!}
\mathfrak{Re}\left[\mathrm{i}^{-n}\,_2F_1\left(\frac{1}{2}, -\frac{n}{2}; \frac{3}{2}; \frac{R^2}{X^2}\right)\right],
\end{split}
}$$

$$\eq{
\begin{split}
\frac{\partial F}{\partial X} = &\phantom{+}\pi e^{-Y} Y_1(X) + \frac{2Y}{XR} \\
& - \frac{2R}{X} e^{-Y} \sum_{n=0}^{\infty} \frac{X^n}{n!}
\mathfrak{Re}\left[\mathrm{i}^{-n}\,_2F_1\left(\frac{1}{2}, -\frac{n}{2}; \frac{3}{2}; \frac{R^2}{X^2}\right)\right],
\end{split}
}$$

$$\eq{
\begin{split}
\frac{\partial^2 F}{\partial X^2} = &\phantom{+}\pi e^{-Y} \left[Y_0(X) - \frac{Y_1(X)}{X} \right] + 
\frac{2Y}{X^2 R} \left(\frac{Y^2}{R^2} - 2 + Y\right) \\
& - \frac{2R}{X^2} e^{-Y} \sum_{n=1}^{\infty} \frac{n X^n}{n!}
\mathfrak{Re}\left[\mathrm{i}^{-n}\,_2F_1\left(\frac{1}{2}, -\frac{n}{2}; \frac{3}{2}; \frac{R^2}{X^2}\right)\right].
\end{split}
}$$

$_2F_1$ is the Gauss hypergeometric function.

### Series Expansion for $D_3$

$$
C_0 = 1 - e^{-Y}, \quad C_n = Y^{2n} - 2n Y^{2n-1} + 2n(2n-1)C_{n-1},
$$

$$\eq{
\begin{split}
F = &-\pi e^{-Y} [H_0(X)+Y_0(X)] - \frac{2(1-e^{-Y})}{X} \\
& - \frac{2}{X} \sum_{n=1}^{\infty} \frac{(-1)^n (2n-1)!!}{2^n n! X^{2n}} C_n,
\end{split}
}$$

$$\eq{
\begin{split}
\frac{\partial F}{\partial X} = &-2e^{-Y} +\pi e^{-Y} [H_1(X)+Y_1(X)] + \frac{2(1-e^{-Y})}{X^2} \\
& + \frac{2}{X^2} \sum_{n=1}^{\infty} \frac{(-1)^n (2n+1) (2n-1)!!}{2^n n! X^{2n}} C_n,
\end{split}
}$$

$$\eq{
\begin{split}
\frac{\partial^2 F}{\partial X^2} = 
&\phantom{+}\pi e^{-Y} \left[H_0(X)+Y_0(X) - \frac{H_1(X)+Y_1(X)}{X}\right] - \frac{4(1-e^{-Y})}{X^3} \\
& - \frac{2}{X^3} \sum_{n=1}^{\infty} \frac{(-1)^n (2n+1) (2n+2) (2n-1)!!}{2^n n! X^{2n}} C_n.
\end{split}
}$$

$(2n-1)!!$ is the odd double factorial.

### Series Expansion for $D_4$

$$
\begin{split}
B_0 = \frac{1-e^{-Y}}{Y}, \quad B_1 = \left(\frac{1}{Y} - \frac{1}{Y^3}\right)e^{-Y} - 2\left(\frac{1}{Y^2} - \frac{1}{Y^3}\right),\\
B_n = \frac{(-1)^{n+1}}{Y}e^{-Y} + \frac{2n(2n-1)}{Y^2}B_{n-1} + \frac{4n(n-1)}{Y^2}B_{n-2},
\end{split}
$$

$$\eq{
\begin{split}
F = &-\pi e^{-Y} [H_0(X)+Y_0(X)] - \frac{2(1-e^{-Y})}{R} \\
& - \frac{2Y}{R} \sum_{n=1}^{\infty} \frac{(-1)^n (2n-1)!!}{2^n n!} \left(\frac{Y}{R}\right)^{2n} B_n,
\end{split}
}$$

$$\eq{
\begin{split}
\frac{\partial F}{\partial X} = &-2e^{-Y} +\pi e^{-Y} [H_1(X)+Y_1(X)] + \frac{2X(1-e^{-Y})}{R^3} \\
& - \frac{2XY}{R^3} \sum_{n=1}^{\infty} \frac{(-1)^n (2n+1)!!}{2^n n!} \left(\frac{Y}{R}\right)^{2n} B_n,
\end{split}
}$$

$$\eq{
\begin{split}
\frac{\partial^2 F}{\partial X^2} = &\phantom{+}\pi e^{-Y} \left[H_0(X)+Y_0(X) - \frac{H_1(X)+Y_1(X)}{X}\right] \\
& + \frac{2(1-e^{-Y})}{R^3} \left(1 - \frac{3X^2}{R^2}\right) \\
& + \frac{2Y}{R^3} \sum_{n=1}^{\infty} \frac{(-1)^n (2n+1)!!}{2^n n!}
\left(\frac{Y}{R}\right)^{2n} \left(1 - \frac{X^2 (2n+3)}{R^2}\right) B_n.
\end{split}
}$$

## Results
The series expansions defined in the previous topic were implemented in python and evaluated for points $(X,Y) \in [0.1, 40] \times [0.1, 40]$. Special functions were evaluated by the [mpmath] library with 16 digits of decimal precision. The images below display $F$, $F_X$ and $F_{XX}$, respectively.

{{< figure src="images/surff.svg" alt="F Surface Plot" align="center" >}}

{{< figure src="images/surffx.svg" alt="Fx Surface Plot" align="center" >}}

{{< figure src="images/surffxx.svg" alt="Fxx Surface Plot" align="center" >}}

The results obtained from the series expansions were compared to the integral expessions $(2)$, $(3)$ and $(4)$. The Integrals were evaluated with [`scipy.quad`][scipyquad] configured with 10 digits of decimal precision. From the following figures, which show the absolute error of $F$, $F_X$ and $F_{XX}$, it is possible to see that the series expansions guarantee at least 8 digits of decimal precision.

{{< figure src="images/errorf.svg" alt="F Absolute error" align="center" >}}

{{< figure src="images/errorfx.svg" alt="Fx Absolute error" align="center" >}}

{{< figure src="images/errorfxx.svg" alt="Fxx Absolute error" align="center" >}}

The integral expressions $(2)$-$(4)$ are 17 times faster than the series expansions because the integrals are actually calculated by the Fortran QUADPACK library. The series expansions, on the other hand, were all implemented in Python. Therefore, the next step is implementing the series expansions in Fortran for better performance.

## Conclusion
This reproduction is a successful first step in the development of the solution of the floating body problem in three dimensions. The next step is to code the series expansions in Fortran. For that, the required special functions are already implemented in [ColSpecF].

## References
1. Penghao Shan and Jiameng Wu. Highly precise approximation of free surface Green function and its high order derivatives based on refined subdomains. Brodogradnja, vol. 69, no. 1, pp. 53–70, 2018. https://doi.org/10.21278/brod69104

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[scipyquad]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
[mpmath]: https://mpmath.org/
[colspecf]: ../../projects/colspecf/
