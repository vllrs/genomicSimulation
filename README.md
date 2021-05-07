# genomicSimulation
*[GitHub Pages link](https://kiravill.github.io/genomicSimulation/)*

This project is an R package that runs stochastic simulations of genetics in breeding schemes. The core of the package is written in C.

This package is not yet submitted to CRAN. [See the latest release here](https://github.com/KiraVill/genomicSimulation/releases).

### Documentation and Guides
[The package vignette is available in the repository at this link](https://kiravill.github.io/genomicSimulation/doc/gSvignette.html). The vignette provides usage guides and an overview of the functions and simulation methodology of the package. 

The docs folder contains the documentation files. These can also be viewed in R using the command `?[function name]`. 

### Additional help files
[A guide to installing package (specifically for the UQ HPC system) is available here.](vignettes/gSinstallguide.md)

### C version
[The C functions that form the base of this package are also provided on GitHub under the name genomicSimulationC.](https://github.com/KiraVill/genomicSimulationC) This version allows you to run the simulation using C function calls.