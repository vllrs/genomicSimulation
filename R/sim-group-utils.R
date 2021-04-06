#' Delete a group's worth of genotypes
#'
#' \code{delete.group} finds the genotypes of the given group and removes
#' them from the SimData object's storage. A message is printed explaining
#' how many genotypes were deleted.
#'
#' @param group an integer representing the group number of the group to be deleted
#' @return 0 on success. An error is raised on failure.
#'
#' @family grouping functions
#' @useDynLib genomicSimulation SXP_delete_group
#' @export
delete.group <- function(group) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_delete_group, sim.data$p, group))
}

#' Assign multiple groups' worth of genotypes to a single group
#'
#' \code{combine.group} combines the members of multple groups
#' into a single group. This is useful for allowing certain crossing
#' functions to act on these genotypes together.
#'
#' @param groups A vector containing the group numbers of groups 
#' that are to be combined
#' @return the group number of the combined group
#'
#' @family grouping functions
#' @useDynLib genomicSimulation SXP_combine_groups
#' @export
combine.groups <- function(groups) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_combine_groups, sim.data$p, length(groups), groups))
}

#' Assign each genotype inside a group to a separate new group
#'
#' \code{break.group.into.individuals} allocates every genotype inside
#' a group to a separate group.
#'
#' @param group an integer: the group number of the group to be split
#' @return the group numbers of the new groups
#'
#' @family grouping functions
#' @useDynLib genomicSimulation SXP_split_individuals
#' @export
break.group.into.individuals <- function(group) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_split_individuals, sim.data$p, group))
}

#' Assign each family inside a group to separate new groups
#'
#' \code{break.group.into.families} allocates all genotypes in the
#' group that share the same parents to separate new groups. 
#'
#' @param group an integer: the group number of the group to be split
#' @return the group numbers of the new groups
#'
#' @family grouping functions
#' @useDynLib genomicSimulation SXP_split_familywise
#' @export
break.group.into.families <- function(group) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_split_familywise, sim.data$p, group))
}

#' Assign certain genotypes to a new group
#'
#' \code{make.group} allocates the genotypes with indexes in the vector
#' indexes to a single new group.
#'
#' @param indexes A vector of integers containing the indexes of the genotypes
#' to move to the new group. Use \code{get.group.data} to get these
#' @return the group number of the new group
#'
#' @family grouping functions
#' @useDynLib genomicSimulation SXP_split_out
#' @export
make.group <- function(indexes) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_split_out, sim.data$p, length(indexes), indexes))
}

#' Get a list of the groups currently existing in the SimData
#'
#' \code{see.existing.groups} collects all groups that currently have genotypes
#' allocated to them and returns the group numbers as a vector.
#'
#' @return A vector containing the numbers of every group in the SimData that 
#' currently has members
#'
#' @family grouping functions
#' @useDynLib genomicSimulation SXP_get_groups
#' @export
see.existing.groups <- function() {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	d <- data.frame(.Call(SXP_get_groups, sim.data$p))
	colnames(d) <- c("Group", "Group size")
	return(d)
}


#' Get data of all individuals that belong to a group.
#'
#' \code{see.group.data} allows you to extract the names or ids or indexes
#' or genotypes of group members as an R vector. Printing to a file then 
#' loading the data from the tab-separated text file is still the suggested
#' method for accessing data, but this may be used to peek at the data or 
#' make group adjustments during simulation time. 
#'
#' When loading genotypes as a group, some genotypes may have trailing nonsense
#' characters. Only the first (number of SNPs * 2) characters are valid. These trailing
#' characters occur because the C code does not allocate extra room for a null-terminating
#' byte, which makes the R function converting the genotypes to strings not know where the
#' end of the genotype is.
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
#' @family saving functions
#' @useDynLib genomicSimulation SXP_get_group_data
#' @export
see.group.data <- function(group, data.type) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_get_group_data, sim.data$p, group, data.type))
}

#' Perform selection on a group by true GEBV.
#'
#' \code{select.by.gebv} calculates the GEBVs of all members of a group, sorts them
#' in either ascending or descending order, and then selects the best few. The number
#' selected can be an exact number (select the best 5) or a percentage (take the best
#' 2\%). Exactly one of percentage and number parameters must be defined for the
#' function to run. 
#'
#' @section Custom selection:
#' For more complex selection than by true/unmasked GEBV, use the function
#' \code{\link{see.group.gebvs}}, use R scripting to identify which to select,
#' and split them into a new group using \code{\link{make.group}}. 
#'
#' @param from.group an integer: the group number of the group perform selection on.
#' @param low.score.best If FALSE, then the highest GEBVs are selected. If TRUE,
#' then the lower ones are selected.
#' @param percentage If this is a number and number and threshold are NULL, then 
#' the [percentage]\% with the best GEBVs will be selected into the new group. Note that
#' this is a whole-number percentage, i.e. for 5\%, set percentage=5 not 0.05.
#' @param number If this is an integer and percentage and threshold are NULL, then
#' the [number] with the best GEBVs will be selected into the new group.
#' @return the group number of the new group, containing the genotypes that were 
#' positively selected.
#'
#' @family grouping functions
#' @useDynLib genomicSimulation SXP_simple_selection
#' @useDynLib genomicSimulation SXP_simple_selection_bypercent
#' @export
select.by.gebv <- function(from.group, low.score.best=FALSE, percentage=NULL, number=NULL) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	if (sum(is.null(percentage), is.null(number)) == 1) {
		if (is.null(percentage)) {
			# we are selecting a certain number
			return(.Call(SXP_simple_selection, sim.data$p, length(from.group), from.group, 
						number, low.score.best))
		} else {
			# we are selecting the top percentage
			return(.Call(SXP_simple_selection_bypercent, sim.data$p, length(from.group), 
						from.group, percentage, low.score.best))
		}
	}
	stop("Exactly one of parameters `percentage` and `number` must be set.")
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
#' @returns a dataframe whose columns are the GEBVs of the group members and 
#' the indexes of the group members.
#'
#' @family grouping functions
#' @family saving functions
#' @useDynLib genomicSimulation SXP_group_eval
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

#' Get a string containing the ultimate/highest-scoring genotype.
#'
#' \code{see.optimal.genotype} allows you to extract the optimal genotype 
#' according to the current loaded effect values as an R string. An error is 
#' raised if no effect values have been loaded.
#'
#' @return A string containing the allele out of available alleles at that
#' marker that has the highest effect value. The string will be ordered in
#' genome order (lowest chromosome and lowest position to highest) according
#' to the map that was included on initialisation.
#'
#' @family saving functions
#' @useDynLib genomicSimulation SXP_get_best_genotype
#' @export
see.optimal.genotype <- function() {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_get_best_genotype, sim.data$p))
}

