% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim-progression.R
\name{self.n.times}
\alias{self.n.times}
\title{Performs n selfing steps on the lines in a group}
\usage{
self.n.times(
  group,
  n,
  map = 0L,
  offspring = 1,
  retain = TRUE,
  give.names = FALSE,
  name.prefix = NULL,
  track.pedigree = TRUE,
  give.ids = TRUE,
  file.prefix = "",
  save.pedigree = FALSE,
  save.gebv = FALSE,
  save.genotype = FALSE
)
}
\arguments{
\item{group}{The group number from which to draw the parents of these
random crosses. If a vector, the parameters and breeding method 
(eg. random crossing) will be applied to each group number in the vector.}

\item{n}{An integer representing the number of steps of selfing to perform.}

\item{map}{The identifier for the recombination map with which gametes 
from members of @a group will be generated. By default uses the oldest 
loaded map currently active in simulation. If a vector of groups is provided,
this parameter can be a single map to be used for crossing operations on all groups, 
or this parameter can be a vector of the same length as `group`, in which case 
each map will be used for simulating the crossing operations on the group which
is at the corresponding position in the vector.}

\item{offspring}{This many offspring of each group member will be produced at the 
first selfing step. After that step, exactly one selfed offspring from each
parent will be progressed, no matter this value.}

\item{retain}{A logical value, repesenting whether to save the generated
genotypes to memory or discard them. You may wish to discard them
but save to file if you are generating too many crosses to save into
memory.}

\item{give.names}{A logical value representing whether or not to produce names
for the new genotypes generated. The names produced would have format [name.prefix][id]}

\item{name.prefix}{A string. If give.names is TRUE, the id is concatenated to 
this to produce the name of each new genotype.}

\item{track.pedigree}{A logical value representing whether or not to save the ids
of the parents of each new genotype to the new genotype's pedigree. If this is 
false, the new genotype's pedigree is unknown.}

\item{give.ids}{A logical value representing whether or not to allocate each new 
genotype an id. If this is FALSE, the new genotype is 'invisible' to pedigree trackers and
even if the pedigree of its offspring is supposedly tracked, the pedigree trackers
will not be able to identify the progenitors of its offspring. Furthermore, if it is
false and names are generated using give.names, all names generated in the same group will
be the same. Probably you'd only have this FALSE if you were discarding the results or worried
about id overflow.}

\item{file.prefix}{A string representing the prefix of files produced if save.pedigree=TRUE, 
save.gebv>0, or save.genotype=TRUE. If no file.prefix is provided, but at least one of 
those settings is TRUE, the default file prefix is "out".
save.pedigree, save.gebv, and save.genotype are known
as the "save-as-you-go" settings, which allow you to output data on the genotypes at the same 
time as they are being generated. The save-as-you-go settings can be used in conjunction with
retain=FALSE, if you want to generate results from a greater number of offspring than your 
system memory can hold.}

\item{save.pedigree}{A logical value. If TRUE, saves the pedigrees of each
generated genotype to the file with filename "[file.prefix]-pedigree.txt". 
Pedigrees are saved progressively (1000 at a time). They are saved in recursive 
predigree format - see @a save.pedigrees or the package vignette for more details.}

\item{save.gebv}{An integer representing a marker effect set identifier, or 0.
If FALSE, 0, or negative, does not save GEBVs. If positive and corresponding to the identifier
of a loaded set of marker effects, saves the GEBVs of each generated genotype to the
file with filename "[file.prefix]-bv.txt". This is a tab-separated text file. 
Breeding values for generated genotypes are calculated and saved progressively (1000 at a 
time). See @a save.GEBVs for more details on the output format. A reminder that 
if only one set of marker effects have been loaded, the identifier for that set of 
marker effects is 1.}

\item{save.genotype}{A logical value. If TRUE, saves the SNP matrix of the generated 
genotypes to the file with filename "[file.prefix]-genotype.txt". 
Generated genotypes are saved progressively (1000 at a time). The output file 
uses genetic markers as columns and genotypes as rows. See @a save.genotypes}
}
\value{
The group number of the new genotypes produced, or 0 if none could be
produced due to an invalid parent group number being provided.
}
\description{
Crosses each genotype with itself repeatedly for a given number of generations.
Returns the group number of the group
that the selfed progeny genotypes were loaded into.
}
\details{
The function performs a selfing step (crosses each genotype
with itself), then performs a selfing step on those offspring, and repeats for a total of 
n generations of selfing. 

The offspring parameter represents the number of genotypes produced from each 
genotype in the original group. These genotypes from the same family are 
independently descended from the original genotype.
}
\seealso{
Other crossing functions: 
\code{\link{make.all.unidirectional.crosses}()},
\code{\link{make.clones}()},
\code{\link{make.crosses.from.file}()},
\code{\link{make.double.crosses.from.file}()},
\code{\link{make.doubled.haploids}()},
\code{\link{make.random.crosses.between}()},
\code{\link{make.random.crosses}()},
\code{\link{make.targeted.crosses}()}
}
\concept{crossing functions}
