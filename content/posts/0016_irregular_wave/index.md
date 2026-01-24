---
date: '2026-01-23'
draft: false
title: 'Generating random waves from a sea spectrum'
author: 'Rodrigo Castro'
summary: 'A long-crested irregular sea is generated from the ISSC spectral formulation.'
tags: ['Water Waves']
---

## Introduction
This post presents a simple method of generating random waves from a sea spectrum. The guideline follows what is presented in the book of *Faltinsen* (1993).

## Methods
Random waves can be seen as a linear combination of waves of different frequencies, whith their amplitudes being regulated by a sea spectrum. This is better explained in the next subtopics.

### Irregular sea
The wave elevation $\zeta$ of a long-crested irregular sea propagating along the positive x-axis can be written as a sum of a large number $N$ of wave components, that is

$$\eq{
\zeta = \sum_{j=1}^{N} A_j \sin (\omega_j t - k_j x + \epsilon_j),
}$$

The frequency $\omega_j$ and the wavenumber $k_j$ are related by the dispersion relation, which for deep waters is expressed as $k_j = \omega_j^2/g$, where $g$ is the acceleration of gravity. The phase angle $\epsilon_j$ is a random number between $0$ and $2\pi$. The amplitude $A_j$ is calculated from the sea spectrum $S(\omega)$ as follows:

$$\eq{
A_j = \sqrt{2 S(\omega_j) \Delta \omega},
}$$

where $\Delta \omega$ is calculated as 

$$\eq{
\Delta \omega = \frac{\omega_{max} - \omega_{min}}{N},
}$$

and $(\omega_{min}, \omega_{max})$ define the range of frequencies that compose the sea. The frequency components $\omega_j$ are chosen randomly in the intervals $[\omega_i, \omega_i + \Delta \omega)$, such that

$$\eq{
\omega_i = \omega_{min} + (i-1) \Delta \omega, \quad i = 1,\ldots,N
}$$

### ISSC spectrum
The wave spectrum describes the energy density as a function of the frequency, and it relies on the assumption that the sea is a stationary random process. The <abbr title="International Ship and Offshore Structures Congress">ISSC</abbr> spectral formulation is one of the most used spectra and is given as

$$\eq{
\frac{S(\omega)}{H_s^2 T_1} = 
\frac{0.11}{2\pi} \left(\frac{\omega T_1}{2\pi}\right)^{-5} \mathrm{exp} \left[-0.44 \left(\frac{\omega T_1}{2\pi}\right)^{-4}\right],
}$$

where $H_s$ is the significant wave height and $T_1$ is the mean wave period. The mean wave period $T_2$ and the modal period $T_0$ are related to $T_1$ by

$$\eq{
T_1 = 1.086 T_2, \quad T_0 = 1.408 T_2.
}$$

## Results
The next image displays the dimensionless ISSC wave spectrum, for $H_s = 8$ and $T_2 = 10$.

{{< figure src="images/issc_h8_t10.svg" alt="ISSC Spectrum" align="center" >}}

With the spectrum above, considering the circular frequency to be in the range $(0.2, 3.2)$, which corresponds to a range of periods between 2.0 and 30 seconds, and 1000 wave components, the following realization (time series) is generated.

{{< figure src="images/ts_h8_t10.svg" alt="Irregular Wave" align="center" >}}

The picture above displays the wave elevation at a point $x=0$. The next animation presents another realization for $x \in [0, 50]$ and a time duration of 60 seconds.

{{< figure src="images/irregular_wave_cropped.gif" alt="Irregular Wave Animation" align="center" >}}

To check that the random wave generator actually follows the starting ISSC spectrum, the Python library [Oceanlyz] is used to do the reverse job. For that, an irregular wave with a time duration of 3 hours and sample rate of 2 Hz is created and given as input. The next image shows the original spectrum against the one obtained by Oceanlyz, with little difference between the two. Also, the significant wave height obtained by Oceanlyz is $H_s=7.9$ and the mean wave period is $T_2=10.1$, very close to the values that were given to generate the ISSC spectrum.

{{< figure src="images/wave_spectrum.svg" alt="Wave Spectrum" align="center" >}}

## Conclusion
The random wave generator works just fine, it can easily be expanded to include other sea spectra, and might be a valuable tool for other articles. One future post will certanly be about the reverse job, that is, how to obtain the spectrum from the wave data.

## References
1. O. M. Faltinsen. 1993. Sea loads on ships and offshore structures. Cambridge University Press, Cambridge, UK.
2. Arash Karimpour and Qin Chen. 2017. Wind wave analysis in depth limited water using OCEANLYZ, a MATLAB toolbox. Computers & Geosciences 106 (Sept. 2017), 181–189. https://doi.org/10.1016/j.cageo.2017.06.010

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[oceanlyz]: https://github.com/akarimp/Oceanlyz
