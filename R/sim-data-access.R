#' Get data of all individuals that belong to a group.
#'
#' \code{see.group.data} allows you to extract the names or ids or indexes
#' or genotypes of group members as an R vector. Printing to a file then 
#' loading the data from the tab-separated text file is still the suggested
#' method for accessing data, but this may be used to peek at the data or 
#' make group adjustments during simulation time. 
#'
#' @section Warning about indexes:
#' Deleting a group (\code{\link{delete.group}}) may cause the indexes of group members to 
#' change (as the positions/indexes that they are stored in the internal matrix may be
#' rearranged). Therefore, you must ensure that there is no call to \code{\link{delete.group}}
#' between loading indexes with this function and passing them to a function that takes indexes
#' as inputs (eg \code{\link{make.group}}, \code{\link{cross}}). 
#'
#' @param group an integer: the group number of the group to have a look at
#' @param data.type a string that will be used to identify the type of data
#' to be extracted. If the first character is 'N' or 'n',
#' the 'N'ames of group members will be extracted. If the first character is 
#' 'D' or 'd', the I'D's of group members will be extracted. If the first character is
#' 'I' or 'i', the 'I'ndexes of group members will be extracted. If the first character
#' is 'G' or 'g', the 'G'enotypes of the group members will be extracted.
#' @return A vector containing the desired data output for each member of the group.
#'
#' @family grouping functions
#' @family data access functions
#' @export
see.group.data <- function(group, data.type) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_get_group_data, sim.data$p, group, data.type))
}


#' Get GEBVs of a group
#'
#' \code{see.group.gebvs} gives the R user a dataframe containing the GEBVs and a matching
#' index for each genotype in a group. This can be used in conjunction with \code{\link{make.group}}
#' to perform a custom selection on a group.
#'
#' @section Usual workflow:
#' The suggested workflow involves using this function to get GEBVs and matching indexes,
#' performing some transformation on the GEBVs (eg environmental masking then
#' sorting), identifying some of the entries for selection, and passing
#' the indexes of those entries to \code{\link{make.group}} in order to split
#' them off into a new group.
#'
#' @param group an integer: the group number of the group to return the GEBVs of.
#' @returns a dataframe whose columns are the GEBVs of the group members, named "GEBV", and 
#' the indexes of the group members, named "i".
#'
#' @family grouping functions
#' @family data access functions
#' @export
see.group.gebvs <- function(group) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	evals <- list()
	for (g in group) {
		d <- data.frame(.Call(SXP_group_eval, sim.data$p, group))
		colnames(d) <- c("i","GEBV")
	} 
	return(d)
}


#' Get a string containing the ultimate/highest-scoring set of alleles.
#'
#' \code{see.optimal.haplotype} allows you to extract the optimal haplotype 
#' according to the current loaded effect values as an R string. An error is 
#' raised if no effect values have been loaded. The optimal genotype (assuming
#' only additive allele effects) is just the doubled version of this haplotype.
#'
#' @return A string containing the allele out of available alleles at that
#' marker that has the highest effect value. The string will be ordered in
#' genome order (lowest chromosome and lowest position to highest) according
#' to the map that was included on initialisation.
#'
#' @family data access functions
#' @export
see.optimal.haplotype <- function() {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_get_best_genotype, sim.data$p))
}


#' Get the ultimate/highest-possible GEBV given the current loaded values.
#'
#' \code{see.optimal.GEBV} allows you to extract the optimal GEBV
#' according to the current loaded effect values. An error is 
#' raised if no effect values have been loaded.
#'
#' @return The highest possible GEBV
#'
#' @family data access functions
#' @export
see.optimal.GEBV <- function() {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_get_best_GEBV, sim.data$p))
}

