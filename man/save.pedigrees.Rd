% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim-printers.R
\name{save.pedigrees}
\alias{save.pedigrees}
\title{Save pedigrees from the simulation to a file}
\usage{
save.pedigrees(filename, group = NULL, type = NULL, recursive.format = TRUE)
}
\arguments{
\item{filename}{A string containing a filename to which the output will
be written}

\item{group}{If not set/set to NULL, will print all genotypes.
Otherwise, if a group of that number exists, save only lines that belong
to that group.}

\item{type}{The printing format. It is now preferred you use `recursive.format`
instead of `type`. For `type`, use a string starting with 'R' or 'r' 
to save in recursive format. Use a string
starting with 'P' or 'p' to save in parents-only format.}

\item{recursive.format}{A logical representing the format of the output pedigree.
If TRUE, print the recursive pedigree that contains all ancestors that genomicSimulation 
can link to a candidate in a nested-bracket formatted file. If FALSE, 
print the immediate parents of the candidate in a three-column file.}
}
\value{
0 on success. On failure an error will be raised.
}
\description{
Saves a record of the ancestors of genotypes
for which this data was tracked.
}
\details{
Two formats are available for printing pedigrees. 

The recursive format (recursive.format=TRUE) recursively traces back and prints
all known parents in the ancestry of the genotype. It first prints the
ID of the genotype, followed by a tab, then recursively prints the pedigree
with format [name of genotype, or ID if no name exists]=([parent1],[parent2]).
(each of the [parent1] and [parent2] will follow a similar format if they
are not founder genotypes or other genotypes with unknown/untracked pedigrees.
One row in the produced file corresponds to the pedigree of one line. 
Genotypes produced by selfing or doubling haploids only show one parent 
(eg. F2=(F1), not F2=(F1,F1)). This is the default printing format.

Parents-only format (recursive.format=FALSE) prints the line's name (or ID if it has none), 
then the two parents' names (or IDs if they have none) to a tab-separated file.
Even if both parents are the same this format prints two columns of parent names. 
Each row of the produced  file is the pedigree of a different genotype. 

Note that the parents-only format prints only the most immediate generation of 
parents of a line. Also note that these 'immediate parents'
are the genotypes from which the ones to be printed were generated (i.e. 
the immediate parents of the result of a call to \code{\link{self.n.times}} are 
the group passed to the \code{\link{self.n.times}} function, not the (n-1)th 
generation of selfing).

The file produced by this function when type="P" can be used as the 
parental.file parameter to \code{\link{find.crossovers}}.
}
\seealso{
Other saving functions: 
\code{\link{change.allele.symbol}()},
\code{\link{save.GEBVs}()},
\code{\link{save.allele.counts}()},
\code{\link{save.genome.model}()},
\code{\link{save.genotypes}()},
\code{\link{save.local.GEBVs.blocks.from.chrsplit}()},
\code{\link{save.local.GEBVs.blocks.from.file}()}
}
\concept{saving functions}
