# genomicSimulation
*[GitHub Pages link](https://vllrs.github.io/genomicSimulation/)*

This project is an R package that runs stochastic simulations of genetics in breeding schemes. The core of the package is written in C. [An overview of the tool, in journal paper format, is available at this link.](http://dx.doi.org/10.1093/g3journal/jkac216).



### Installation

This package is not yet available on CRAN. 

[See the latest release here](https://github.com/vllrs/genomicSimulation/releases). Alternatively, install the development version from [Github](https://github.com/vllrs/genomicSimulation).

```r
remotes::install_github('vllrs/genomicSimulation')
```


### Documentation and Guides

An introduction and guides to simulating with the tool are available in [the package vignette](https://vllrs.github.io/genomicSimulation/doc/gSvignette.html). Alternatively, the  worked examples in the [Features and Template Guide (for genomicSimulation and genomicSimulationC)](https://vllrs.github.io/genomicSimulationC/html/templates.html) can be used as a different introduction to the package.  

The docs folder of this repository contains the documentation files for every function in the package. This documentation can also be viewed in R using the command `?[function name]`. 

### Additional help files
[A guide to installing package (specifically for the UQ HPC system) is available here.](doc/gSinstallguide.md)

### C version
[The C functions that form the base of this package are also provided on GitHub under the name genomicSimulationC.](https://github.com/vllrs/genomicSimulationC) This version allows you to run the simulation using C function calls.
