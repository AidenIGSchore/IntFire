### Aiden I. G. Schore
### Department of Plant Biology, University of Illinois at Urbana-Champaign
### [aschore2@illinois.edu](mailto:aschore2@illinois.edu)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14360231.svg)](https://doi.org/10.5281/zenodo.14360231)

This repository contains the code used in "Plant Functional types improve satellite-derived burn severity assessments in interior Alaska" to calculate
the new Integrated Vegetation-Fire Severity Index (IntFire).

The R code compares satellite based burn indices with ground based burn severity measurements using linear and logarithmic regressions.
The spreadsheets with the data used with the R code in the paper can be found at the [Arctic Data Center](https://doi.org/10.18739/A25H7BW3S).

'Fire Unmixing' is a copy of the [FireUnmixing](https://code.earthengine.google.com/460e8e7e862bf6f3fc72c8ec3e585786) Google Earth Engine script for the project. 
Using the link is recommended over cloning the provided file.

After running the first two sections of the R code and the Google Earth Engine script, the final sections of the R code run validation of IntFire's performance.

Reference to be added upon publication of paper.
