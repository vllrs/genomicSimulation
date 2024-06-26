% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim-group-utils.R
\name{break.group.into.buckets}
\alias{break.group.into.buckets}
\title{Randomly split all genotypes in a group between a set of specified-capacity
buckets}
\usage{
break.group.into.buckets(group, buckets)
}
\arguments{
\item{group}{an integer: the group number of the group to be split}

\item{buckets}{a vector of integers: the number of genotypes to be allocated 
to the first, second, third, etc. buckets.}
}
\value{
the group numbers of length(buckets)+1 groups into which the genotypes were split
}
\description{
Given a `buckets` parameter of length n, n+1 groups will be created, with the 
last group containing (number of group members) - (sum of other bucket 
capacities) genotypes. That is, it is assumed that the sum of the input 
capacities is less than the number of group members, and the leftover amount 
goes in the last bucket.
}
\details{
If the bucket capacities add up to more than the number of group members, 
a warning will be raised, but the function will still run. Note, though, that
buckets are filled first to last. 
If the total capacity of the set of buckets is greater than the number of 
group members, later buckets in the list will not be filled to capacity.
}
\seealso{
Other grouping functions: 
\code{\link{break.group.by.GEBV}()},
\code{\link{break.group.by.label.range}()},
\code{\link{break.group.by.label.value}()},
\code{\link{break.group.evenly}()},
\code{\link{break.group.into.families}()},
\code{\link{break.group.into.halfsib.families}()},
\code{\link{break.group.into.individuals}()},
\code{\link{break.group.randomly}()},
\code{\link{break.group.with.probabilities}()},
\code{\link{change.label.default}()},
\code{\link{combine.groups}()},
\code{\link{delete.group}()},
\code{\link{make.group}()},
\code{\link{see.existing.groups}()},
\code{\link{see.group.data}()},
\code{\link{see.group.gene.data}()}
}
\concept{grouping functions}
