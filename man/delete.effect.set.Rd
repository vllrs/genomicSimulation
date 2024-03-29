% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim-deletors.R
\name{delete.effect.set}
\alias{delete.effect.set}
\title{Delete marker effect set(s)}
\usage{
delete.effect.set(effect.sets)
}
\arguments{
\item{effect.sets}{an vector containing the label numbers of the labels to be deleted}
}
\value{
0 on success. An error is raised on failure.
}
\description{
\code{delete.effect.set} removes the given marker effect set(s) from simulation memory.
}
\seealso{
Other grouping functions: 
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
\code{\link{delete.label}()},
\code{\link{make.group.from.label.range}()},
\code{\link{make.group.from.label}()},
\code{\link{make.group}()},
\code{\link{see.existing.groups}()},
\code{\link{see.group.data}()},
\code{\link{see.group.gene.data}()},
\code{\link{select.by.gebv}()}
}
\concept{grouping functions}
