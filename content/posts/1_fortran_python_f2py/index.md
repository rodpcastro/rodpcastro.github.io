---
date: '2025-04-26'
draft: false
title: 'Integrating Fortran and Python with F2PY'
author: 'Rodrigo Castro'
summary: 'Calling Fortran procedures from Python.'
tags: ['Fortran', 'Python', 'F2PY']
---

## Introduction
With future projects demanding mathematical [special functions], I’ve [started][spfuncs] implementing them in [Fortran] for its computational efficiency. However, not everything will be built with pure Fortran, as I also plan to integrate it with the versatility of Python through F2PY. So the purpose of this post is to demonstrate how to bridge Fortran and Python with F2PY.

## Methods

### Tools
[F2PY] is a tool included with the [NumPy] package. To function, it requires a Fortran compiler, such as [gfortran]. Furthermore, for Python versions 3.12 and above, F2PY demands the installation of [meson] backend due to [changes in NumPy's build process][NumPy migration].

### Fortran code
The Fortran subroutine to be wrapped by F2PY contains the numerical implementation of the exponential integral:

$$ \text{Ei}(x) = \int_{-\infty}^x \frac{e^t}{t} \, dt, \quad x \in \mathbb{R}, x > 0 $$

The expressions for numerical evaluation of \\(\text{Ei}\\) and many other special functions can be found on the book *[Computation of Special Functions][csf_book]* by Zhang and Jin. The numerical implementation given below is a refit of the code found in the book.

<details><summary>expint.f90</summary>
{{< include file="expint.md" >}}
</details>

The `expi` subroutine might not display the most correct way of writing subroutines in Fortran, but it reflects the adaptations suggested by [F2PY examples] for correct wrapping:

- *Avoid Fortran modules*: In a typical Fortran setup, we might place the `expi` subroutine inside a module for organization. However, since F2PY automatically wraps the subroutine into a Python module, this would create a redundant nested structure. To keep the Python interface clean, the subroutine is defined outside of a Fortran module.

- *Handling `iso_fortran_env` parameters*: `int16` and `real64` had to be used inside the subroutine. Unclear why, these parameters were the cause of parsing errors in F2PY when used outside of the subroutine (within a module).

### F2PY
To ensure that F2PY correctly interprets Fortran `int16` and `real64` kinds, they need to be [mapped][f2cmap] to their equivalent C types. This is done by using a `.f2py_c2map` file placed in the same directory as the `expint.f90` file:

```python
# .f2py_c2map
dict(real=dict(real64='double'), integer=dict(int16='int'))
```

After that, all that remains is to build our extension module by running the following command:

```console
f2py -c expint.f90 -m expint
```
<details><summary><i>shell output</i></summary>
{{< include file="out.md" >}}
</details>

The command produces a shared library (`expint.cp312-win_amd64.pyd`). The `expi` procedure can now be called from a Python module named `expint`.

## Results
To demonstrate that the `expi` procedure works as expected, its outputs are compared against [mpmath], a Python library for arbitrary precision arithmetic. The `mpmath.ei` function is used with 16 digits of precision and the relative error between `expi` and `mpmath.ei` is computed and displayed through the following Python code:

<details><summary><i>expint_test.py</i></summary>
{{< include file="expint_test.md" >}}
</details>

The resulting plot shows that the relative difference between `expi` and `mpmath.ei` is very small (≤ \\(10^{-13}\\)) for \\(10^{-10} ≤ x ≤ 10^{2}\\).

{{< figure src="expi_error.svg" alt="expi relative error" align="center" >}}

## References
1. NumPy Developers. 2024. F2PY. Distributed as part of NumPy. https://numpy.org/doc/stable/f2py/index.html. (2025).

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[Fortran]: https://fortran-lang.org/
[NumPy]: https://numpy.org/
[F2PY]: https://numpy.org/doc/stable/f2py/
[gfortran]: https://gcc.gnu.org/wiki/GFortran
[meson]: https://github.com/mesonbuild/meson
[mpmath]: https://mpmath.org/
[special functions]: https://en.wikipedia.org/wiki/Special_functions
[exponential integral]: https://en.wikipedia.org/wiki/Exponential_integral
[spfuncs]: https://github.com/rodpcastro/special-functions
[csf_book]: https://search.worldcat.org/title/33971114
[F2PY examples]: https://numpy.org/doc/stable/f2py/f2py-examples.html
[f2cmap]: https://numpy.org/doc/stable/f2py/advanced/use_cases.html#dealing-with-kind-specifiers
[NumPy migration]: https://numpy.org/doc/stable/reference/distutils_status_migration.html#distutils-status-migration
