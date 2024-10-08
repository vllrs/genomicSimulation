% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim-printers.R
\name{save.genome.model}
\alias{save.genome.model}
\title{Save the details of the genome that the SimData uses.}
\usage{
save.genome.model(filename)
}
\arguments{
\item{filename}{A string containing a filename to which the output will
be written}
}
\value{
0 on success. On failure an error will be raised.
}
\description{
**DEFUNCT**
> This function was never correctly implemented. 
An implementation may return at some future date.
Original documentation below.
}
\details{
Saves the SNP names, their linkage map positions, 
and the effect values for calculating GEBVs (if applicable) to a file.

The output is printed with first a tab-separated header row, containing
the column names "name", "chr", "pos", then all alleles for which the saved
SimData has effect values, if it has any. The subsequent lines, also tab-separated,
contain the name of a SNP, the chromosome on which it is found, the positions
at which it is found, and then the contribution to the GEBV of every allele
at that marker for which effects are loaded.
}
\seealso{
Other saving functions: 
\code{\link{change.allele.symbol}()},
\code{\link{save.GEBVs}()},
\code{\link{save.allele.counts}()},
\code{\link{save.genotypes}()},
\code{\link{save.local.GEBVs.blocks.from.chrsplit}()},
\code{\link{save.local.GEBVs.blocks.from.file}()},
\code{\link{save.pedigrees}()}
}
\concept{saving functions}
\keyword{internal}
