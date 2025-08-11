---
date: '2025-06-16'
draft: false
title: 'Accelerating Python code with Numba and F2PY'
author: 'Rodrigo Castro'
summary: 'Using Numba and F2PY to speed up computations in Python.'
tags: ['Python', 'Numba', 'Fortran', 'F2PY']
---

## Introduction
Following the suggestion mentioned in the conclusion of the [previous article][2d_incoef], the objective of this post is to speed up the computation of the 2D constant boundary element influence coefficients, $\mathbb{G}$ and $\mathbb{Q}$, by two methods and compare their performances.

## Methods
The two methods to accelerate Python code that are studied in this post are [Numba] and [F2PY]. They are further examined in the following sections.

### Numba
[Numba] is a open-source just-in-time compiler that translates a subset of Python and NumPy into fast machine code. The most common way of using Numba is through its decorators, which can be applied to the functions the user wants to compile. The following code snippet displays the application of Numba decorator [`jit`][numba.jit] to a class static method that computes the influence coefficients $\mathbb{G}$ and $\mathbb{Q}$ by analytical integration.

{{< dropdown_file title="Class static method accelerated with Numba" src="analytical_numba.py" fmt="python" >}}

The decorator could have been applied to a function, but it was used to a static method of the class `Element` since this calculation is exclusive to its objects. In Numba version 0.61.0, applying the decorator to an instance method was not possible, as the `jit` compiler could not interpret the `self` argument.

The used arguments to `numba.jit` are described below: 

* `'Tuple((f8, f8))(f8[:], f8)'`: this is the signature that represents the types of the arguments, where `f8` means 8-byte float number. The first part `Tuple((f8, f8))` indicates two float outputs, the second part `(f8[:], f8)` describe one float array and one float variable as inputs. In general, the signature is optional, but it's necessary when `cache=True` is used.
* `nopython=True`: this argument tells Numba to run the static method entirely without the involvement of the Python interpreter. This is the recommended practice for best performance.
* `cache=True`: this option makes Numba use a cached version of the compiled method. This avoids recompilation every time the routine is called for the first time.

### F2PY
[F2PY], already explored in [another post][post_f2py], is a tool that builds Python wrappers for Fortran routines, allowing the speed of Fortran compiled code to be used from within the Python environment. Below is the Fortran code to be wrapped by F2PY.

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
    
    # [...]
        
    def get_influence_coefficients_gauss_fortran(self, field_global):
        field_local = self.get_point_local_coordinates(field_global)

        return incoef.incoef.gauss_fortran(field_local, self.length)
```

In the code above, the first `incoef` is the Python module created by F2PY. The second `incoef` is the Fortran module wrapped by F2PY. The inclusion of the Fortran subroutine inside a Fortran module allows the addition of auxilary functions that can be kept private from the Python environment.

## Results
The [previous post][2d_incoef] compared analytical integration and the 4-point Gauss-Legendre quadrature. Here, each method is independently reassessed, enhanced by the previously mentioned Numba and F2PY techniques, with focus on their performance improvements. The following plot shows the computation time for calculating the influence coefficients $\mathbb{G}$ and $\mathbb{Q}$ for several field points layed around a boundary element, the same conditions presented in the previous post.

{{< figure src="computation_time.svg" alt="Computation time" align="center">}}

On average, methods improved by Numba and Fortran perform similary, taking <ins>half the time of the original procedures</ins>. 

## Conclusion
The results indicate that both Numba and Fortran significantly enhance a Python routine performance. Numba's ease of use encourages its further application, though it's best suited for numerically oriented code using basic Python elements and NumPy, with limitations on other Python constructs. Conversely, Fortran offers greater flexibility in optimizing Python code, but the F2PY wrapper imposes constraints, requiring Fortran modules and procedures to follow a strict format distinct from traditional Fortran development.

## References
1. Siu Kwan Lam, Antoine Pitrou, and Stanley Seibert. 2015. Numba: a LLVM-based Python JIT compiler. In Proceedings of the Second Workshop on the LLVM Compiler Infrastructure in HPC (LLVM '15). Association for Computing Machinery, New York, NY, USA, Article 7, 1–6. https://doi.org/10.1145/2833157.2833162
2. NumPy Developers. 2024. F2PY. Distributed as part of NumPy. https://numpy.org/doc/stable/f2py/

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[2d_incoef]: ../0003_2d_constant_boundary_element/
[numba]: https://numba.pydata.org/
[f2py]: https://numpy.org/doc/stable/f2py/
[numba.jit]: https://numba.pydata.org/numba-doc/dev/reference/jit-compilation.html#numba.jit
[post_f2py]: ../0001_f2py_fortran_python/
