---
date: '2026-06-05'
draft: true
title: 'Intact stability analysis of a fishing vessel'
author: 'Rodrigo Castro'
summary: 'The intact stability of a fishing vessel, in an unique loading condition, is studied and checked against the International Code on Intact Stability.'
tags: ['Intact Stability', 'navaltoolbox', 'FreeCAD']
---

## Introduction
In this post is presented a basic intact stability analysis of a fishing vessel, which modeling was described in the [last article][fctrawler]. The vessel geometry is built with [FreeCAD] and the stability calculations are performed by [NavalToolbox][nvtbox]. The results are checked against the [International Code on Intact Stability][iscode].

## Methods
The analysis starts at the CAD modeling stage, to correctly determine what geometrical features can contribute to the vessel's stability. The next step is estimating the vessel's mass properties: displacement and center of gravity. Then, the geometry and mass properties are used by a specialized software to get basic hydrostatic properties and perform the actual intact stability analysis, which results are validated (or not) by criteria defined in international rules.

### 3D modeling
The surface that was modeled in the [previous post][fctrawler] has to be modified and incremented in order to get all the geometrical features that can contribute to the ship's buoyancy. First, that surface is cut along the main deck level. Next, cabins, forecastle and wheelhouse are added. The final elements that make up all the enclosed volumes are presented in the following image.

{{< figure src="images/trawler_closed_volumes.png" alt="Closed volumes" align="center" >}}

The doors that give access to the cabins and wheelhouse are weathertight, therefore, their sills must be identified as downflooding points. The next image shows these points' positions.

{{< figure src="images/trawler_weathertight_doors.png" alt="Weathertight doors" align="center" >}}

Another important aspect of the ship is its windage, which consists of the projected area of the vessel that is exposed to the wind. This area is represented by the red silhouette drawn in the next figure. The points that define this red boundary are exported in a `csv` file to be used by the stability software.

{{< figure src="images/trawler_windage.png" alt="Windage" align="center" >}}

The geometry that represents the enclosed volumes of the vessel must be given to the stability software as a `stl` mesh. Thankfully, FreeCAD has [Gmsh] and other meshing tools integrated, which easily allowed the creation of the following mesh, containing 3602 faces.

{{< figure src="images/trawler_mesh.png" alt="Mesh" align="center" >}}

### Estimation of mass properties
The ship is studied at its design water level, therefore, displacement and longitudinal position of the center of gravity (LCG) are already known, or can be easily obtained with [NavalToolbox][nvtbox]'s [`HydrostaticCalculator`][nvtbox_fdraft]. No initial trim or list is assumed, so LCG = LCB (longitudinal center of buoyancy), and the transversal center of gravity (TCG) is zero. The remaining of the task is estimaing the vertical position of the center of gravity (KG). For this, we use the work of *Kim & Yeo* (2020). In this reference, the KG can be calculated as a function of the gross tonnage, for multiple loading conditions typical of fishing vessels. In the present article, it is assumed that the design waterline coincides with what *Kim & Yeo* define as "Full load departure", which configures the condition where the ship is fully loaded with consumables such as fuel, fresh water and provisions.

There is one caveat to this estimation process. When reading the article of *Kim & Yeo*, it's noticeable, and verifiable by simple calculations, that the volumes occupied by cabins and wheelhouses are not included in the gross tonnage of the fishing vessels used as reference in their work. For this reason, the present work considers as enclosed volume, for the purpose of gross tonnage calculation, only the volume of the hull below the main deck level.

### Intact stability analysis


## Results


## Conclusion


## References
1. Food and Agriculture Organization. 2026. Fishing Vessel Design Database - Trawler - 17.5m. https://www.fao.org/fishery/en/vesseldesign/mar-17
2. Dong Jin Kim, Dong Jin Yeo. 2020. Estimation of drafts and metacentric heights of small fishing vessels according to loading conditions. International Journal of Naval Architecture and Ocean Engineering 12 (2020), 199–212. https://doi.org/10.1016/j.ijnaoe.2019.11.001
3. Antoine Anceau. 2026. NavalToolbox Version 0.9.0. https://github.com/NavalToolbox/navaltoolbox-lib.
4. International Maritime Organization. 2008. International Code on Intact Stability. International Maritime Organization, London.

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[fctrawler]: ../0019_3dcad_trawler-mor_freecad
[freecad]: https://www.freecad.org/
[gmsh]: https://gmsh.info/
[nvtbox]: https://navaltoolbox.github.io/navaltoolbox-lib/
[nvtbox_fdraft]: https://navaltoolbox.github.io/navaltoolbox-lib/api/python/hydrostatics.html#HydrostaticsCalculator.from_draft
[iscode]: https://wwwcdn.imo.org/localresources/en/KnowledgeCentre/IndexofIMOResolutions/MSCResolutions/MSC.267(85).pdf
