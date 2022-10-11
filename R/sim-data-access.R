#' Get data of all individuals that belong to a group.
#'
#' \code{see.group.data} allows you to extract the names, IDs, indexes, 
#' genotypes (as strings), breeding values, or parentage of group members as an R vector.
#' As long as the group is not modified in between calls, the same ordering of group
#' members will be used, so this function can be used to construct tables of data, eg.
#' \code{data.frame("i"=see.group.data(3,"X"),"GEBV"=see.group.data(3,"BV"))}.
#'
#' @section Warning about indexes:
#' Deleting a group (\code{\link{delete.group}}) may cause the indexes of group members to 
#' change (as the positions/indexes that they are stored in the internal matrix may be
#' rearranged). Therefore, you must ensure that there is no call to \code{\link{delete.group}}
#' between loading indexes with this function and passing them to a function that takes indexes
#' as inputs (eg \code{\link{make.group}}).
#'
#' @param group an integer: the group number of the group to have a look at
#' @param data.type a string that will be used to identify the type of data
#' to be extracted. String parsing is case-insensitive.
#' If the first character is 'N', the 'N'ames of group members will be extracted. 
#' If the first character is 'D', the I'D's of group members will be extracted. 
#' If the first character is 'X', the Inde'X'es of group members will be extracted. 
#' If the first character is 'G', the 'G'enotypes of the group members will be extracted. 
#' If the first character is 'B', the 'b'reeding values/GEBVs of group members will be returned.
#' If it starts with "P1", the first parent of each group member is returned. As usual in 
#' genomicSimulation, if the parent has a name, its name will be what is returned, otherwise,
#' it will be its ID.
#' If it starts with "P2", the second parent of each group member is returned. As usual in 
#' genomicSimulation, if the parent has a name, its name will be what is returned, otherwise,
#' it will be its ID.
#' If it starts with "ped", the full pedigree log of each group member is returned.
#' @return A vector containing the desired data output for each member of the group.
#'
#' @family grouping functions
#' @family data access functions
#' @export
see.group.data <- function(group, data.type) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_get_group_data, sim.data$p, group, toupper(data.type)))
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
	return(.Call(SXP_get_best_haplotype, sim.data$p))
}

#' Get a string containing the ultimate/highest-scoring set of alleles available
#' in the current group.
#' 
#' That is, consider the pool of alleles that exist in the group for each locus,
#' and take the highest-scoring allele at each locus.
#'
#' An error is raised if no effect values have been loaded.
#'
#' @param group an integer: the group number of the group to have a look at. 
#' Can be a vector of groups, in which case the calculation will be run independently
#' for each group.
#' @return A string containing the allele out of all alleles in the group
#' that has the highest effect value at that locus. The string will be ordered in
#' genome order (lowest chromosome and lowest position to highest) according
#' to the map that was included on initialisation.
#'
#' @family data access functions
#' @export
see.optimal.possible.haplotype <- function(group) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_get_best_available_haplotype, sim.data$p, length(group), group))  
}

#' Get the ultimate/highest-possible GEBV given the current loaded effect values.
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

#' Get the ultimate/highest-possible GEBV given the pool of alleles available in
#' the current group.
#' 
#' Effectively, given the additive model of trait effects, this is 2 times the 
#' breeding value score of \code{see.optimal.possible.haplotype}. 
#'
#' An error is raised if no effect values have been loaded.
#'
#' @param group an integer: the group number of the group to have a look at.
#' Can be a vector of groups, in which case the calculation will be run independently
#' for each group.
#' @return The highest possible GEBV that can be created from alleles available
#' in the group.
#'
#' @family data access functions
#' @export
see.optimal.possible.GEBV <- function(group) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_get_best_available_GEBV, sim.data$p,length(group), group))    
}

#' Get the lowest-possible GEBV given the current loaded effect values.
#'
#' \code{see.minimum.GEBV} allows you to extract the lowest GEBV score
#' possible according to the current loaded effect values. An error is 
#' raised if no effect values have been loaded.
#'
#' @return The lowest possible GEBV
#'
#' @family data access functions
#' @export
see.minimum.GEBV <- function() {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_get_worst_GEBV, sim.data$p))
}