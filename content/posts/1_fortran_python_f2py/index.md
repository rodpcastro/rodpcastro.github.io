---
date: '2025-04-25'
draft: false
title: 'Integrating Fortran and Python with F2PY'
author: 'Rodrigo Castro'
summary: 'Calling Fortran procedures from Python.'
tags: ['Fortran', 'Python', 'F2PY']
---

## Introduction
With future projects demanding mathematical [special functions], I’ve started [implementing][spfuncs] them in [Fortran] for its computational efficiency. However, not everything will be built with pure Fortran – I also plan to integrate it with the versatility of Python through [F2PY]. So the purpose of this post is to demonstrate how to bridge Fortran and Python with F2PY.

## Methods

### Tools
[F2PY] is a tool installed with the [NumPy] package. To have it working, it's also necessary the installation of a Fortran compiler, like [gfortran]. Furthermore, for Python ≥ 3.12, F2PY requires the installation of [meson] backend.

### Fortran code
The Fortran subroutine to be wrapped by F2PY contains the numerical implementation of the exponential integral:

$$ Ei(x) = -\int_{-x}^{\infty} \frac{e^{-t}}{t} dt $$

The expressions for numerical evaluation of \\(Ei\\) and many other special functions can be found on the book [Computation of special functions], by Zhang and Jin. Following the recommendations of this book, the numerical implementation of the exponential integral in Fortran is given below:

<details><summary>expint.f90</summary>
{{< include file="expint.md" >}}
</details>

Reclecting the [F2PY examples], the code had to be adapted for the generation of correct wrapping. In a Fortran setup, the subroutine would have been placed inside of a module. However, if this was adopted here, the `expi` procedure would unecessarily be positioned inside two modules after wrapping. Also, it was noted that, the call for the `iso_fortran_env` parameters `int16` and `real64` raised errors when placed outside of the subroutine, despite correct [mapping between Fortran kinds and C types](#f2py).

### F2PY
To Fortran kinds `int16` and `real64` be interpreted by F2PY, they need be mapped to their equivalent C types. This is done through the file `.f2py_c2map` placed in the same folder as the `.f90` file. Adopting [F2PY definitions][f2cmap], our mapping file contains:

```python
dict(real=dict(real64='double'), integer=dict(int16='int'))
```

After that, all that remains is to build our extension module by running:

```console
f2py -c expint.f90 -m expint
```
<details><summary><i>shell output</i></summary>
{{< include file="out.md" >}}
</details>

The [resulting wrapper](#appendices) (`.pyd`) containing the `expi` procedure can now be called from Python.

## Result
To demonstrate that our new Fortran procedure works as expected, the results of our exponential integral implementation are going to be compared against [mpmath], a solid Python library that can compute the value of the exponential integral with arbitrary precision.

## Appendices
* <a href="expint.f90" download>expint.f90</a>
* <a href="expint.cp312-win_amd64.pyd" download>expint.cp312-win_amd64.pyd</a>
* <a href="expint_test.py" download>expint_test.py</a>

<!--Links-->
[Fortran]: https://www.manning.com/books/modern-fortran
[NumPy]: https://numpy.org/
[F2PY]: https://numpy.org/doc/stable/f2py/
[gfortran]: https://gcc.gnu.org/wiki/GFortran
[meson]: https://github.com/mesonbuild/meson
[mpmath]: https://mpmath.org/
[special functions]: https://en.wikipedia.org/wiki/Special_functions
[exponential integral]: https://en.wikipedia.org/wiki/Exponential_integral
[spfuncs]: https://github.com/rodpcastro/special-functions
[Computation of special functions]: https://search.worldcat.org/title/33971114
[F2PY examples]: https://numpy.org/doc/stable/f2py/f2py-examples.html
[f2cmap]: https://numpy.org/doc/stable/f2py/advanced/use_cases.html
