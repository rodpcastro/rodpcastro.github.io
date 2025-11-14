---
date: '2025-11-14'
draft: false
title: 'GDeepSeries'
author: 'Rodrigo Castro'
summary: 'Fortran implementation of the 3D infinite-depth free-surface Green function according to series expansions.'
tags: ['Green Function', 'Fortran', 'Potential Flow']
---

## Introduction
This Fortran library computes the three-dimensional infinite-depth free-surface Green function according to the series expansions defined in the work of [*Shan* and *Wu* (2018)](#reference).

An overview of the expressions implemented in this library can be found [here][src]. Alternatively, a Python implementation of the series expansions can be found in this [blog post][rpcgds].

## Tests
Simple tests were performed to compare the values computed by GDeepSeries against the numerical evaluation with Python libraries. The test results are summarized [here][test].

## Dependency
GDeepSeries makes use of special functions implemented in [ColSpecF].

## Reference
1. Penghao Shan and Jiameng Wu. Highly precise approximation of free surface Green function and its high order derivatives based on refined subdomains. Brodogradnja, vol. 69, no. 1, pp. 53–70, 2018. <https://doi.org/10.21278/brod69104>

## License
The GDeepSeries code is distributed under the MIT License (see [LICENSE] file).

**Important Dependency Notice**: This project depends on [ColSpecF], which is distributed under multiple licenses. Review [ColSpecF's license][csf-license] for details.

<!-- links -->
[colspecf]: ../colspecf/
[src]: https://github.com/rodpcastro/gdeepseries/tree/main/src/README.md#infinite-depth-free-surface-green-function
[license]: https://github.com/rodpcastro/gdeepseries/blob/main/LICENSE
[csf-license]: https://github.com/rodpcastro/colspecf/blob/main/LICENSE
[rpcgds]: ../../posts/0013_3d_inf_depth_fsurface_gfunction/
[test]: https://github.com/rodpcastro/gdeepseries/blob/main/test/README.md#test-results
