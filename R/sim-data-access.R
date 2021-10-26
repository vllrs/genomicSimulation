#' Get data of all individuals that belong to a group.
#'
#' \code{see.group.data} allows you to extract the names, IDs, indexes, 
#' genotypes (as strings), breeding values, or parentage of group members as an R vector.
#' As long as the group is not modified in between calls, the same ordering of group
#' members will be used, so this function can be used to construct tables of data, eg.
#' \code{data.frame("i"=see.group.data(3,"X"),"GEBV"=see.group.data(3,"BV"))}.
#' 
#' @Section Deprecation:
#' To avoid confusion, the data.type parameter code for indexes has been changed from 'I' to
#' 'X' to reduce ambiguity with the data type IDs ('D'). 
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
  if (substr(toupper(data.type),1,1) == "I") {
    warning("Data type code \"I\" is deprecated; to get group indexes, please use \"X\".", call. = FALSE)
    return(.Call(SXP_get_group_data, sim.data$p, group, "X"))
  }
	return(.Call(SXP_get_group_data, sim.data$p, group, toupper(data.type)))
}


#' Get GEBVs of a group
#'
#' \code{see.group.gebvs} gives the R user a dataframe containing the GEBVs and a matching
#' index for each genotype in a group. This can be used in conjunction with \code{\link{make.group}}
#' to perform a custom selection on a group.
#'
#' @section Deprecation:
#' In the history of this package, there was no way to get GEBVs using \code{\link{see.group.data}},
#' because that function only dealt with data that required no calculation step. However, it 
#' now includes that capability, and this function has become nothing but a wrapper. Deprecating it
#' and exposing the underlying call:
#' \code{data.frame("i"=see.group.data(...,data.type="X"),"GEBV"=see.group.data(...,data.type="BV")))}
#' demystifies this function and allows users to see how they can use \code{\link{see.group.data}} for 
#' flexible data access.
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
  .Deprecated("data.frame(\"i\"=see.group.data(...,\"X\"),\"GEBV\"=see.group.data(...,\"BV\"))")
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
	return(.Call(SXP_get_best_haplotype, sim.data$p))
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