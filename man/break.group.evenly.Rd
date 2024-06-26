% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim-group-utils.R
\name{break.group.evenly}
\alias{break.group.evenly}
\title{Randomly split all genotypes in a group between n identical-size groups.}
\usage{
break.group.evenly(group, into.n = 2)
}
\arguments{
\item{group}{an integer: the group number of the group to be split}

\item{into.n}{an integer: the number of groups into which to break this one}
}
\value{
the group numbers of the new groups
}
\description{
The sizes of the resulting groups will differ by at most 1, in the case 
where n does not divide the number of group members perfectly.
}
\seealso{
Other grouping functions: 
\code{\link{break.group.by.GEBV}()},
\code{\link{break.group.by.label.range}()},
\code{\link{break.group.by.label.value}()},
\code{\link{break.group.into.buckets}()},
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
