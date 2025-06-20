---
date: '2025-06-09'
lastmod: '2025-06-20'
draft: false
title: 'ColSpecF'
author: 'Rodrigo Castro'
summary: 'ColSpecF is a Fortran library for evaluating mathematical Special Functions.'
tags: ['Special Functions', 'Fortran']
---

{{< figure src="colspecf_logo.svg" alt="ColSpecF Logo" align="center">}}

## Introduction
[ColSpecF][ColSpecF GitHub] (Collected Special Functions) is a [Fortran][Fortran Website] library for evaluating mathematical [Special Functions], built around adaptations of algorithms collected from [several sources](#references).

## Tests
Tests are conducted by comparing the ColSpecF results with those of [mpmath], an arbitrary-precision numerical library.

Testing routines are built using [test-drive], a standard Fortran unit testing framework.

## Documentation
The [API documentation][ColSpecF Docs] for this library is generated using [FORD] and is deployed and hosted on [ReadTheDocs].

## References
Fortran code for evaluating special functions is sourced from the following websites:

* Association for Computing Machinery. 2012. [Collected Algorithms][calgo]
* Jason Blevins. 2004. [Alan Miller's Fortran Software][jblevins]
* Commonwealth Scientific and Industrial Research Organisation. 2004. [Software from Alan J. Miller][csiro]
* Elsevier. 2025. [Elsevier Data Repository][elsvdata]
* John Burkardt. 2025. [Fortran77 Source Codes][jbf77]
* John Burkardt. 2025. [Fortran90 Codes][jbf90]

This list grows as more special functions are added to the library.

## License
ColSpecF is a Fortran library distributed under multiple licenses or permissions based on code origin. Users must comply with the applicable license or permission for each portion of the code. See the [License][License File] for full details.

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
[calgo]: https://calgo.acm.org/
[jblevins]: https://jblevins.org/mirror/amiller/
[csiro]: https://wp.csiro.au/alanmiller/
[elsvdata]: https://elsevier.digitalcommonsdata.com/
[jbf77]: https://people.sc.fsu.edu/~jburkardt/f77_src/f77_src.html
[jbf90]: https://people.sc.fsu.edu/~jburkardt/f_src/f_src.html
<!-- License -->
[License File]: https://github.com/rodpcastro/colspecf/blob/main/LICENSE
