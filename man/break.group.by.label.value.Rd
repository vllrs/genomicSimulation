% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim-group-utils.R
\name{break.group.by.label.value}
\alias{break.group.by.label.value}
\title{Assign genotypes with a certain label value to a new group}
\usage{
break.group.by.label.value(label, value, group = NA)
}
\arguments{
\item{label}{the label number of the label to look at.}

\item{value}{all genotypes selected by the `group` parameter that have 
this integer as their label value of the appropriate label will be moved to
the new group.}

\item{group}{NA or 0 to search through all genotypes in simulation memory,
or a group number or vector of group numbers to search only
members of those groups.}
}
\value{
the group number of the new group
}
\description{
Allocates the genotypes with a particular value
for a given custom label to a new group.
}
\details{
Multiple group inputs supported. All selected genotypes, no matter their 
group of origin, will be moved into the same output group.
}
\seealso{
Other label functions: 
\code{\link{break.group.by.label.range}()},
\code{\link{change.label.by.amount}()},
\code{\link{change.label.to.this}()},
\code{\link{change.label.to.values}()},
\code{\link{change.names.to.values}()},
\code{\link{create.new.label}()},
\code{\link{delete.label}()}

Other grouping functions: 
\code{\link{break.group.by.GEBV}()},
\code{\link{break.group.by.label.range}()},
\code{\link{break.group.evenly}()},
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
\concept{label functions}
