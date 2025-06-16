---
date: '2025-06-15'
draft: true
title: 'Accelerating Python code with Numba and F2PY'
author: 'Rodrigo Castro'
summary: 'Using Numba and F2PY to speed up computations in Python.'
tags: ['Python', 'Numba', 'Fortran', 'F2PY']
---

## Introduction
Following the suggestion mentioned in the conclusion of the [previous article][2d_incoef], the objective of this post is to speed up the computation of the 2D constant boundary element influence coefficients, $\mathbb{G}$ and $\mathbb{Q}$, by two methods and compare their performances.

## Methods
The two approaches to accelerate Python code that are studied in this post are [Numba] and [F2PY]. They are further explained in the subsequent topics.

### Numba
[Numba] is a open-source just-in-time compiler that translates a subset of Python and NumPy into fast machine code. The most common way of using Numba is through its decorators, which can be applied to the functions the user wants to compile. The following code snippet displays the application of Numba decorator [`jit`][numba.jit] to a class static method that computes the influence coefficients $\mathbb{G}$ and $\mathbb{Q}$ by analytical integration.

{{< dropdown_file title="Class static method accelerated with Numba" src="analytical_numba.py" fmt="python" >}}

The same decorator could have been applied to a function, but the choice to decorate a static method of the class `Element` is because this calculation is only performed for objects of this class. In the used version of Numba (`0.61.0`), the same decorator could not be applied to an instance method, because the argument `self` could not be interpreted by the `jit` compiler.

In the code above the arguments to `numba.jit` represent:

* `'Tuple((f8, f8))(f8[:], f8)'`: this is the signature that represents the types of the method's arguments, where `f8` means 8-byte float number. The string indicates two float outputs `Tuple((f8, f8))`, one float array and one float number as inputs (`f8[:], f8`). In general, this input is optional, but it's necessary when `cache=True` is used.
* `nopython=True`: this argument tells Numba to run the static method entirely without the involvement of the Python interpreter. This is the recommended and best-practice way to use Numba `jit` decorator as it leads to the best performance.
* `cache=True`: this option makes Numba use a cached version of the compiled method. This avoids recompilation every time the routine is called for the first time.

### F2PY
[F2PY] was already explored in [this post][post_f2py]. F2PY is a tool that builds Python wrappers for Fortran routines, allowing the speed of Fortran compiled code to be used from within the Python environment. Below is the Fortran code to be wrapped by F2PY.

{{< dropdown_file title="Fortran subroutine" src="gauss_fortran.f90" fmt="fortran" >}}

The `incoef` python module is created by running:

```console
f2py -c incoef.f90 -m incoef
```

After creation of a shared library (`incoef.cp312-win_amd64.pyd` on the author machine), the `gauss_fortran` subroutine can then be called from a Python instance method like this:

```python
# Extract from element.ipynb

import incoef

class Element:
    
    [...]
        
    def get_influence_coefficients_gauss_fortran(self, field_global):
        field_local = self.get_point_local_coordinates(field_global)

        return incoef.incoef.gauss_fortran(field_local, self.length)
```

In the code above, the first `incoef` is the Python module created by F2PY. The second `incoef` is the Fortran module wrapped by F2PY. The inclusion of the Fortran subroutine inside a Fortran module allows the addition of auxilary functions that can be kept private from the Python environment.

## Results
In the [previous post][2d_incoef], analytical integration and the 4-point Gauss-Legendre quadrature were compared in terms of performance. Here they are reevaluated separately, each improved by the two methods aforementioned, and now the focus is on how much the two original procedures are improved by Numba and Fortran (F2PY). The following plot displays the computation time taken to calculate the influence coefficients $\mathbb{G}$ and $\mathbb{Q}$ for several field points layed around a boundary element, the same conditions presented in the [previous post][2d_incoef].

{{< figure src="computation_time.svg" alt="Computation time" align="center">}}

On average, both Numba and Fortran perform similary, taking half the time of the original procedures. 


## Conclusion


## References
1. 
2. NumPy Developers. 2024. F2PY. Distributed as part of NumPy. https://numpy.org/doc/stable/f2py/index.html. (2025).

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[2d_incoef]: ../3_2d_constant_boundary_element/
[numba]: https://numba.pydata.org/
[f2py]: https://numpy.org/doc/stable/f2py/
[numba.jit]: https://numba.pydata.org/numba-doc/dev/reference/jit-compilation.html#numba.jit
[post_f2py]: ../1_f2py_fortran_python/
