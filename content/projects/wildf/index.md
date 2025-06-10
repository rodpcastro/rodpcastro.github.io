---
date: '2025-06-09'
draft: false 
title: 'WildF'
author: 'Rodrigo Castro'
summary: 'WildF is a Fortran library for evaluating mathematical Special Functions.'
tags: ['Special Functions', 'Fortran']
---

{{< figure src="wildf_logo.svg" alt="WildF Logo" align="center">}}

## Introduction
[WildF][WildF GitHub] is a [Fortran][Fortran Website] library for evaluating mathematical [Special Functions]. Just like wild creatures, Special Functions are unusual, but interesting mathematical objects that arise in many areas of applied mathematics. This library aims to serve as a support for any project that needs the computation of these rare species.

## Tests
Tests are conducted by comparing the WildF results with those of [mpmath], an arbitrary-precision numerical library. These tests ensure at least 8 digits of precision within the specified domains.

Testing routines are built using [test-drive], a standard Fortran unit testing framework.

## Documentation
The [API documentation][WildF Docs] for this library is generated using [FORD] and is deployed and hosted on [ReadTheDocs].

## References
1. Shanjie Zhang, Jianming Jin. 1996. [Computation of Special Functions][Book Zhang]. Wiley, New York, NY.

<!-- Links -->
<!-- Introduction -->
[WildF GitHub]: https://github.com/rodpcastro/wildf
[Fortran Website]: https://fortran-lang.org/
[Special Functions]: https://www.britannica.com/science/special-function
<!-- Tests -->
[mpmath]: https://mpmath.org/
[test-drive]: https://github.com/fortran-lang/test-drive
<!-- Documentation -->
[WildF Docs]: https://wildf.readthedocs.io/
[FORD]: https://forddocs.readthedocs.io/
[ReadTheDocs]: https://about.readthedocs.com/
<!-- References -->
[Book Zhang]: https://search.worldcat.org/title/33971114
