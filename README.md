# genomicSimulation
*[GitHub Pages link](https://vllrs.github.io/genomicSimulation/)*

This project is an R package that runs stochastic simulations of genetics in breeding schemes. The core of the package is written in C. [An overview of the tool, in journal paper format, is available at this link.](http://dx.doi.org/10.1093/g3journal/jkac216).


### Installation

#### Latest Release
This package is not yet available on CRAN. You can [find the latest release here](https://github.com/vllrs/genomicSimulation/releases), and install it using the following command:

```r
install.packages('genomicSimulation_0.2.6.zip', repos=NULL)
# or 
install.packages('genomicSimulation_0.2.6.tar.gz', repos=NULL)
```

The binary (.zip) install is compiled for Windows and R version 4.3, and may not be compatible with all systems. If the binary package cannot be installed on your system, use the source package (.tar.gz). Installing the source package requires C compilers.

If you are from UQ and installing the package on bunya, there is a specific guide for that [at this link.](doc/gSinstallguide-bunya.md)

#### The Very Latest Features
Alternatively, install the development version from [Github](https://github.com/vllrs/genomicSimulation). This also requires C compilers. 

```r
remotes::install_github('vllrs/genomicSimulation')
```

### Documentation and Guides

An introduction and guides to simulating with the tool are available in [the package vignette](https://vllrs.github.io/genomicSimulation/doc/gSvignette.html). Alternatively, the  worked examples in the [Features and Template Guide (for genomicSimulation and genomicSimulationC)](https://vllrs.github.io/genomicSimulationC/html/templates.html) can be used as a different introduction to the package.  

The docs folder of this repository contains the documentation files for every function in the package. This documentation can also be viewed in R using the command `?[function name]`. 

#### Additional help
If you have questions about the package, do open an issue on GitHub or send an email to the package maintainer email address.

### C version
[The C functions that form the base of this package are also provided on GitHub under the name genomicSimulationC.](https://github.com/vllrs/genomicSimulationC)
