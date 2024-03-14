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
	return(.Call(SXP_combine_groups, sim.data$p, groups))
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
	return(.Call(SXP_break_group_into_individuals, sim.data$p, group))
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
	return(.Call(SXP_break_group_into_families, sim.data$p, group))
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
  return(.Call(SXP_break_group_into_halfsib_families, sim.data$p, group, parent))
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
  return(.Call(SXP_break_group_randomly, sim.data$p, group, into.n))  
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
  return(.Call(SXP_break_group_evenly, sim.data$p, group, into.n))  
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
  return(.Call(SXP_break_group_into_buckets, sim.data$p, group, buckets))  
}

#' Randomly split all genotypes in a group between a set of buckets with provided 
#' probabilities.
#' 
#' Given a `probabilities` parameter of length n, n+1 groups will be created, with 
#'  1 - (sum of other probabilities) probability of being in the last group. If 
#' probabilities add up to 1, a warning will be raised but no extra group will 
#' be created.
#' That is, it is assumed that the sum of `probabilities` is less than 1, 
#' with the remaining probability falling towards a last group. The reason for this 
#' is that floating-point numbers have their imperfect accuracies, so we want to 
#' save the computer from having to guess if you meant a last probability with 
#' 0.000001 probability or if you were trying to provide some decimals that added
#' up nicely to 1. 
#'
#' If the bucket capacities add up to more than 1, 
#' a warning will be raised, but  the function will still run. Note, though, that
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
  return(.Call(SXP_break_group_into_probabilities, sim.data$p, group, probabilities))  
}

#' Assign certain genotypes to a new group by index
#'
#' \code{make.group} allocates the genotypes with indexes in the vector
#' indexes to a single new group.
#'
#' Unlike other functions in this package with the prefix \code{make},
#' \code{make.group} does not create any new genotypes, only changes the 
#' group allocation of already-existing genotypes.
#'
#' @param indexes A vector of integers containing the indexes of the genotypes
#' to move to the new group. Use \code{get.group.data} to get these
#' @return the group number of the new group
#'
#' @family grouping functions
#' @export
make.group <- function(indexes) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_make_group_from, sim.data$p, indexes))
}

#' Assign genotypes with a certain label value to a new group
#'
#' \code{break.group.by.label.value} allocates the genotypes with a particular value
#' for a given custom label to a new group.
#' 
#' Multiple group inputs supported. All selected genotypes, no matter their 
#' group of origin, will be moved into the same output group.
#'
#' @param label the label number of the label to look at.
#' @param value all genotypes selected by the `group` parameter that have 
#' this integer as their label value of the appropriate label will be moved to
#' the new group.
#' @param group NA or 0 to search through all genotypes in simulation memory,
#' or a group number or vector of group numbers to search only
#' members of those groups.
#' @return the group number of the new group
#'
#' @family label functions
#' @family grouping functions
#' @export
break.group.by.label.value <- function(label, value, group=NA) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_break_group_by_label_value, sim.data$p, label, value, group))
}

#' OLD NAME | Assign genotypes with a certain label value to a new group
#' 
#' ! This is the old name for \code{break.group.by.label.value}. From genomicSimulation v0.2.5,
#' \code{break.group.by.label.value} is the recommended name over \code{make.group.from.label}. 
#' \code{make.group.from.label} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{break.group.by.label.value}
#' 
#' @keywords internal 
#' @export
make.group.from.label <- function(label, value, group=NA) {
  return(break.group.by.label.value(label,value,group))
}

#' Assign genotypes with label values in a certain range to a new group
#'
#' \code{break.group.by.label.range} allocates the genotypes with label values
#' for a given custom label in a certain range from `rangeLowEnd` to `rangeHighEnd`
#' inclusive to a new group.
#' 
#' Multiple group inputs supported. All selected genotypes, no matter their 
#' group of origin, will be moved into the same output group.
#'
#' @param label the label number of the label to look at.
#' @param rangeLowEnd minimum integer value of label `label` that is to be put
#' in the new group.
#' @param rangeHighEnd maximum integer value of label `label` that is to be put
#' in the new group.
#' @param group NA or 0 to search through all genotypes in simulation memory,
#' or a group number or vector of group numbers to search only
#' members of those groups.
#' @return the group number of the new group
#'
#' @family label functions
#' @family grouping functions
#' @export
break.group.by.label.range <- function(label, rangeLowEnd, rangeHighEnd, group=NA) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_break_group_by_label_range, sim.data$p, label, rangeLowEnd, rangeHighEnd,
               group))
}

#' OLD NAME | Assign genotypes with label values in a certain range to a new group
#' 
#' ! This is the old name for \code{break.group.by.label.range}. From genomicSimulation v0.2.5,
#' \code{break.group.by.label.range} is the recommended name over \code{make.group.from.label.range}. 
#' \code{make.group.from.label.range} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{break.group.by.label.range}
#' 
#' @keywords internal 
#' @export
make.group.from.label.range <- function(label, rangeLowEnd, rangeHighEnd, group=NA) {
  return(break.group.by.label.range(label,rangeLowEnd,rangeHighEnd,group))
}


#' Change the default value of a custom label
#'
#' \code{change.label.default} changes the default (birth) value of a custom
#' label(s) to the given integer value(s).
#' 
#' If `defaults` is shorter than `labels` 
#'
#' @param labels the label(s) of which to change the default values.
#' @param defaults the new default value of the new custom label. All genotypes
#' generated in future will have this value for this label. If there are multiple
#' labels, then `defaults` should be a vector of the same length, with the new
#' defaults for each label in `labels` at corresponding positions.
#' @return 0 on success
#'
#' @family grouping functions
#' @export
change.label.default <- function(labels, defaults) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_change_label_default, sim.data$p, labels, defaults))
}


#' Set values of a custom label
#'
#' \code{change.label.to.values} changes the values of a custom label for all 
#' genotypes or for members of group(s) to a sequence of values.
#' 
#' The function orders the genotypes selected by the `group` parameter by their
#' index, finds the `startIndex`th genotype (counting up from 0), then copies
#' the vector `values` entry-by-entry into the appropriate label of 
#' the genotypes from there on.
#' 
#' If the vector of values provided is longer than the number of genotypes after
#' the `startIndex`th that are selected by the `group` parameter, then trailing 
#' values are ignored. 
#' If the vector of values provided is shorter than the number of genotypes after
#' the `startIndex`th that are selected by the `group` parameter, then trailing 
#' genotypes are left unchanged. 
#'
#' @param label the label number of the label to change.
#' @param values A vector of integer values to copy into the label. 
#' @param group NA or 0 to apply the values to all genotypes in simulation memory,
#' or a group number to apply the values to only members of that group.
#' @param startIndex The position within the list of genotypes identified by `group`
#' and sorted by their index that the first entry of `values` should be copied into.
#' The first position in the list is considered index 1.
#' @return 0 on success
#'
#' @family label functions
#' @export
change.label.to.values <- function(label, values, group=NA, startIndex=1) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_change_label_to_values, sim.data$p, label, 
               values, group, startIndex))    
}


#' Set a flat value for a custom label
#'
#' \code{change.label.to} changes the value of a custom label for all 
#' genotypes or for members of group(s) to a certain value.
#' 
#' Multiple group input supported.
#'
#' @param label the label number of the label to set.
#' @param value an integer to which the label value of all genotypes selected 
#' by the `group` parameter will be changed.
#' @param group NA or 0 to change all genotypes in simulation memory,
#' or a group number or vector of group numbers to change only
#' members of those groups.
#' @return 0 on success
#'
#' @family label functions
#' @export
change.label.to.this <- function(label, value, group=NA) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_change_label_to_this, sim.data$p, label, value, group))  
}


#' Add/subtract from a custom label
#'
#' \code{change.label.by.amount} changes the values of a custom
#' label for all genotypes or for members of group(s) by a given increment amount, 
#' positive or negative.
#' 
#' Multiple group input supported.
#'
#' @param label the label number of the label to change.
#' @param amount the integer change to modify the label value for the chosen
#' genotypes. eg. 1 for incrementing label values by 1, -3 to subtract 3 from
#' all label values.
#' @param group NA or 0 to apply this increment to all genotypes in simulation memory,
#' or a group number or vector of group numbers to apply the incrememnt to only
#' members of those groups.
#' @return 0 on success
#'
#' @family label functions
#' @export
change.label.by.amount <- function(label, amount, group=NA) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_change_label_by_amount, sim.data$p, label, amount, group))
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
#' @param eff.set identifier for the set of marker effects to use to calculate the
#' GEBVs.
#' @return the group number of the new group, containing the genotypes that were 
#' positively selected.
#'
#' @family grouping functions
#' @export
#' @aliases break.group.by.gebv
break.group.by.GEBV <- function(from.group, low.score.best=FALSE, percentage=NULL, number=NULL, eff.set=1L) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	if (sum(is.null(percentage), is.null(number)) == 1) {
		if (is.null(percentage)) {
			# we are selecting a certain number
			return(.Call(SXP_break_group_by_GEBV_num, sim.data$p, from.group, eff.set,
						number, low.score.best))
		} else {
			# we are selecting the top percentage
			return(.Call(SXP_break_group_by_GEBV_percent, sim.data$p, from.group, eff.set,
			             percentage, low.score.best))
		}
	}
	stop("Exactly one of parameters `percentage` and `number` must be set.")
}

#' OLD NAME | Perform selection on a group by true GEBV.
#' 
#' ! This is the old name for \code{break.group.by.GEBV}. From genomicSimulation v0.2.5,
#' \code{break.group.by.GEBV} is the recommended name over \code{select.by.gebv}. 
#' \code{select.by.gebv} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{break.group.by.GEBV}
#' 
#' @keywords internal 
#' @export
select.by.gebv <- function(from.group, low.score.best=FALSE, percentage=NULL, number=NULL, eff.set=1L) {
  return(break.group.by.GEBV(from.group,low.score.best,percentage,number,eff.set))
}