# Guide to installing package on UQ's bunya system

1. If you don't already have a folder for storing personal R packages set up, do this first:
	1. Make a directory for your R packages. 
		- For example: `mkdir ~/Rlibs`. Replace `~/Rlibs` with your choice of folder name.
	2. Tell R where to find the R packages that will be saved here:
		- `echo '.libPaths( c("~/Rlibs", .libPaths()) )' >> ~/.Rprofile`, with `~/Rlibs` modified to be the folder from step 1.
2. To install this package:
	1. Download the genomicSimulation source package (.tar.gz) from the Releases page on Github (https://github.com/vllrs/genomicSimulation/releases), and save that .tar.gz file somewhere in your directories on bunya.
		- For file transfer between your machine and bunya, the RCC recommends using scp or WinSCP. See their guides.
		- In this tutorial, the downloaded package is assumed to be located at `~/genomicSimulation_0.1.tar.gz` on bunya
	2. Get an interactive job on an "epyc3" node: `salloc --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=5G --job-name=TinyInteractive --time=01:00:00 --partition=general --qos=debug --constraint=epyc3 --account=ACCOUNT srun `, where ACCOUNT is replaced with your bunya group account string. 
		- It's important to compile software on bunya on an epyc3 node. Otherwise, if you are later running a job, and bunya allocates the job to an epyc3 node, the job may crash with an "invalid operand" or "illegal instruction" error when it tries to use a package that was compiled on a non-epyc3 node. 
	3. On the interactive job node, open the R version of your choice.
		- For example, run `module load r`, then open the R commandline by typing `R`
	4. In R, on the interactive job node, install the package from the release file: `install.packages("~/genomicSimulation_0.1.tar.gz", repos=NULL)` with the location you saved the package release substituted for `~/genomicSimulation_0.1.tar.gz`
		- Once this process has succeeded, the `~/genomicSimulation_0.1.tar.gz` file can be deleted. The `.Rprofile` file and the contents of the `~/Rlibs` folder should not be deleted. You can run `q()` to exit R, and `exit` to release your interactive job and return to the login node.

Once the package is installed, you can use `library("genomicSimulation")` in any R script. There should be no need to rerun the `install.packages` command, even after `module unload` and `module load`.

