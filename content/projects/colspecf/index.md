---
date: '2025-06-09'
draft: false
title: 'ColSpecF'
author: 'Rodrigo Castro'
summary: 'ColSpecF is a Fortran library for evaluating mathematical Special Functions.'
tags: ['Special Functions', 'Fortran']
---

{{< figure src="colspecf_logo.svg" alt="ColSpecF Logo" align="center">}}

## Introduction
[ColSpecF][ColSpecF GitHub] (Collected Special Functions) is a [Fortran][Fortran Website] library for evaluating mathematical [Special Functions], built around adaptations of [Collected Algorithms][calgo] from [ACM] to modern Fortran.

## Tests
Tests are conducted by comparing the ColSpecF results with those of [mpmath], an arbitrary-precision numerical library. These tests ensure at least 8 digits of precision within the specified domains.

Testing routines are built using [test-drive], a standard Fortran unit testing framework.

## Documentation
The [API documentation][ColSpecF Docs] for this library is generated using [FORD] and is deployed and hosted on [ReadTheDocs].

## References
Fortran code for evaluating special functions is sourced from the following websites:

1. Association for Computing Machinery. 2012. [Collected Algorithms][calgo]
2. Jason Blevins. 2004. [Alan Miller's Fortran Software][jblevins]
3. Commonwealth Scientific and Industrial Research Organisation. 2004. [Software from Alan J. Miller][csiro]

## License
ColSpecF is distributed under two licenses based on code origin:

- Code comprising original contributions by Rodrigo Castro is licensed under the MIT License.
- Code adapted from Collected Algorithms (CALGO), published by the Association for Computing Machinery (ACM), is subject to the [ACM Software License Agreement][acmlic].

Users must comply with the applicable license for each portion of the code. See the [License][License File] for full details.

<!-- links -->
<!-- Introduction -->
[ColSpecF GitHub]: https://github.com/rodpcastro/colspecf
[Fortran Website]: https://fortran-lang.org/
[Special Functions]: https://www.britannica.com/science/special-function
<!-- Tests -->
[mpmath]: https://mpmath.org/
[test-drive]: https://github.com/fortran-lang/test-drive
<!-- Documentation -->
[ColSpecF Docs]:https://colspecf.readthedocs.io/
[FORD]: https://forddocs.readthedocs.io/
[ReadTheDocs]: https://about.readthedocs.com/
<!-- References -->
[acm]: https://www.acm.org/
[calgo]: https://calgo.acm.org/
[jblevins]: https://jblevins.org/mirror/amiller/
[csiro]: https://wp.csiro.au/alanmiller/
<!-- License -->
[acmlic]: https://www.acm.org/publications/policies/software-copyright-notice
[License File]: https://github.com/rodpcastro/colspecf/blob/main/LICENSE
