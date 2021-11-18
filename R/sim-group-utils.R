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

#' Assign all genotypes inside a group that share one parent 
#' to separate new groups
#'
#' \code{break.group.into.halfsib.families} allocates all genotypes in the
#' group that share one same parent (either the same first or the same 
#' second parent) to separate new groups. 
#'
#' @param group an integer: the group number of the group to be split
#' @param parent an integer: 1 to group by the first parent, 2 to group 
#' by the second
#' @return the group numbers of the new groups
#'
#' @family grouping functions
#' @export
break.group.into.halfsib.families <- function(group, parent) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_split_halfsibwise, sim.data$p, group, parent))
}

#' Randomly assign each genotype inside a group to one of n groups, with equal
#' probability.
#'
#' The sizes of the resulting groups are dependent on the random draws.
#'
#' @param group an integer: the group number of the group to be split
#' @param into.n an integer: the number of groups into which to break this one
#' @return the group numbers of the new groups
#'
#' @family grouping functions
#' @export
break.group.randomly <- function(group, into.n=2) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_split_randomly, sim.data$p, group, into.n))  
}

#' Randomly split all genotypes in a group between n identical-size groups.
#'
#' The sizes of the resulting groups will differ by at most 1, in the case 
#' where n does not divide the number of group members perfectly.
#'
#' @param group an integer: the group number of the group to be split
#' @param into.n an integer: the number of groups into which to break this one
#' @return the group numbers of the new groups
#'
#' @family grouping functions
#' @export
break.group.evenly <- function(group, into.n=2) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_split_evenly, sim.data$p, group, into.n))  
}

#' Randomly split all genotypes in a group between a set of specified-capacity
#' buckets
#' 
#' Given a `buckets` parameter of length n, n+1 groups will be created, with the 
#' last group containing (number of group members) - (sum of other bucket 
#' capacities) genotypes. That is, it is assumed that the sum of the input 
#' capacities is less than the number of group members, and the leftover amount 
#' goes in the last bucket.
#'
#' If the bucket capacities add up to more than the number of group members, 
#' a warning will be raised, but the function will still run. Note, though, that
#' buckets are filled first to last. 
#' If the total capacity of the set of buckets is greater than the number of 
#' group members, later buckets in the list will not be filled to capacity. 
#'
#' @param group an integer: the group number of the group to be split
#' @param buckets a vector of integers: the number of genotypes to be allocated 
#' to the first, second, third, etc. buckets.
#' @return the group numbers of length(buckets)+1 groups into which the genotypes were split
#'
#' @family grouping functions
#' @export
break.group.into.buckets <- function(group, buckets) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_split_buckets, sim.data$p, group, buckets))  
}

#' Randomly split all genotypes in a group between a set of buckets with provided 
#' probabilities.
#' 
#' Given a `probabilities` parameter of length n, n+1 groups will be created, with 
#'  1 - (sum of other probabilities) probability of being in the last group. If 
#' probabilities add up to 1, a warning will be raised but no extra group will 
#' be created.
#' That is, it is assumed that the sum of `probabilities` is less than 1, 
#' with the remaining probability falling towards a last group.
#'
#' If the bucket capacities add up to more than 1, 
#' a warning will be raised, but the function will still run. Note, though, that
#' genotypes are allocated to buckets using cumulative probabilities, so groups 
#' with cumulative probabilities of more than 1 will get no members.
#'
#' @param group an integer: the group number of the group to be split
#' @param probabilities a vector of decimals: the probability a genotypes will be 
#' allocated to the first, second, third, etc. buckets.
#' @return the group numbers of the new groups
#'
#' @family grouping functions
#' @export
break.group.with.probabilities <- function(group, probabilities) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_split_probabilities, sim.data$p, group, probabilities))  
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
#' \code{\link{see.group.data}}, use R scripting to identify which to select,
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
