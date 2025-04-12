---
date: '2025-04-12'
draft: false 
title: 'Nodimo'
author: 'Rodrigo Castro'
tags: ['Python', 'Jupyter Notebook', 'Sympy', 'Dimensional Analysis']
summary: 'Nodimo turns a dimensional relationship between quantities into a dimensionless expression.'
---

<p align="center">
<img src="/images/nodimo.svg" alt="Nodimo">
</p>

# Nodimo
The main purpose of Nodimo is to transform a dimensional relationship between quantities into a dimensionless one. This is done by grouping dimensional quantities into dimensionless products in such a way that the resulting number of products is always lower than or equal to the starting number of quantities. Therefore, the ensuing dimensionless model is, at the same time, a generalization and simplification of the dimensional model.

Nodimo supports any number of dimensions and quantities. It can be used for applications in science, engineering, economics and finance. The resulting dimensionless relations can be used as the basis for further studies in similarity and model testing.

## Notes

* The use of Nodimo requires basic knowledge of dimensional analysis, specially on choosing the appropriate set of scaling parameters and indentifying established dimensionless groups.

* It is recommended the use of [jupyter notebook][Jupyter Notebook] for a better displaying of the results.

## Installation
Via `PyPI`, Nodimo is installed by:
```shell
pip install nodimo
```

Alternatively, via `Anaconda`:
```shell
conda install nodimo
```

## Getting started
### Basic example
* Simple pendulum

<p align="center">
    <img src="/images/simple_pendulum.svg" alt="Simple Pendulum">
</p>

The dimensionless relation between the pendulum's period `T` and the other quantities presented in the figure above is built and displayed as:

```python
from nodimo import Quantity, Model

T = Quantity('T', mass=0, length=0, time=1, dependent=True)  # period
L = Quantity('L', mass=0, length=1, time=0, scaling=True)    # length
m = Quantity('m', mass=1, length=0, time=0)                  # mass
g = Quantity('g', mass=0, length=1, time=-2, scaling=True)   # gravity
t0 = Quantity('theta_0')                                     # initial angle

model = Model(T, L, m, g, t0)
model.show()
```

And the result is:

$$\frac{T g^{\frac{1}{2}}}{L^{\frac{1}{2}}} = \Phi{\left(\theta_{0} \right)}$$

For more applications and functionalities, check the [documentation][Docs Status].

<!-- Links -->
[Docs Status]: https://nodimo.readthedocs.io/
[Jupyter Notebook]: https://jupyter.org/
