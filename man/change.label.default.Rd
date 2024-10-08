% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim-group-utils.R
\name{change.label.default}
\alias{change.label.default}
\title{Change the default value of a custom label}
\usage{
change.label.default(labels, defaults)
}
\arguments{
\item{labels}{the label(s) of which to change the default values.}

\item{defaults}{the new default value of the new custom label. All genotypes
generated in future will have this value for this label. If there are multiple
labels, then `defaults` should be a vector of the same length, with the new
defaults for each label in `labels` at corresponding positions.}
}
\value{
0 on success
}
\description{
Changes the default (birth) value of a custom
label(s) to the given integer value(s).
}
\details{
If `defaults` is shorter than `labels`
}
\seealso{
Other grouping functions: 
\code{\link{break.group.by.GEBV}()},
\code{\link{break.group.by.label.range}()},
\code{\link{break.group.by.label.value}()},
\code{\link{break.group.evenly}()},
\code{\link{break.group.into.buckets}()},
\code{\link{break.group.into.families}()},
\code{\link{break.group.into.halfsib.families}()},
\code{\link{break.group.into.individuals}()},
\code{\link{break.group.randomly}()},
\code{\link{break.group.with.probabilities}()},
\code{\link{combine.groups}()},
\code{\link{delete.group}()},
\code{\link{make.group}()},
\code{\link{see.existing.groups}()},
\code{\link{see.group.data}()},
\code{\link{see.group.gene.data}()}
}
\concept{grouping functions}
