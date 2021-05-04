# Guide to installing package on UQ HPC

1. If you do not have a folder for storing your personal R packages already set up, perform these tasks first:
	1. Make a directory to which your personal R packages can be loaded. 
		- For example: `mkdir ~/Rlibs`
	2. Tell R where to find the R packages that will be saved here:
		- `echo '.libPaths( c("~/Rlibs", .libPaths()) )' >> ~/.Rprofile`, with `~/Rlibs` modified to be the folder from step 1.
2. Steps specific to installing this package:
	1. Download the genomicSimulation release package (.tar.gz) to the HPC system. 
		- The RCC site provides examples of software that can be used to transfer files. Some examples are WinSCP, scp, and FileZilla.
		- For the purposes of the following steps, suppose the package release is saved at `~/genomicSimulation_0.1-1.tar.gz`
	2. Open the R version of your choice.
		- For example, load with `module load R/3.5.0`, then open the R commandline using `R`
	3. Install the package from the release file:
		- `install.packages("~/genomicSimulation_0.1-1.tar.gz", repos=NULL)` with the location you saved the package release substituted for `~/genomicSimulation_0.1-1.tar.gz`
		- Once this process has succeeded, the `~/genomicSimulation_0.1-1.tar.gz` file can be deleted. The `.Rprofile` file and the contents of the `Rlibs` folder should not be deleted.

From this point, the R command `library("genomicSimulation")` will allow you to use the package functions. There should be no need to rerun the `install.packages` command, even after `module unload` and `module load`.
