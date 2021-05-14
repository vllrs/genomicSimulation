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
#' @export
make.group <- function(indexes) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_split_out, sim.data$p, length(indexes), indexes))
}

#' Get a list of the groups currently existing in the SimData
#'
#' \code{see.existing.groups} scans the saved data for groups that currently
#' have members and returns their group numbers and number of members.
#'
#' @return A dataframe containing two columns, the first, named "Group", 
#' being the group numbers of every group in the SimData that currently 
#' has members, the second, named "GroupSize", being 
#' the number of genotypes currently allocated to that group
#'
#' @family grouping functions
#' @export
see.existing.groups <- function() {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	d <- data.frame(.Call(SXP_get_groups, sim.data$p))
	colnames(d) <- c("Group", "GroupSize")
	return(d)
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
