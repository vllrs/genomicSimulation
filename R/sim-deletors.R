#' Clear the internal storage of all data
#'
#' Frees and deletes all data stored in the
#' package's internal SimData struct (sim.data$p). Returns 0 on success.
#'
#' @family loader functions
#' @family deletor functions
#' @export
clear.simdata <- function() {
	if (is.null(sim.data$p)) { return(0L) }
	sim.data$p <- NULL
	gc()
}


#' Delete a group's worth of genotypes
#'
#' Finds the genotypes of the given group and removes
#' them from the SimData object's storage. A message is printed explaining
#' how many genotypes were deleted.
#' 
#' Multiple group input supported.
#'
#' @param groups an vector containing the group numbers of the groups to be deleted
#'
#' @family grouping functions
#' @family deletor functions
#' @export
delete.group <- function(groups) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  groups <- convert.to.integer(groups,"Group identifiers")
	.Call(SXP_delete_group, sim.data$p, groups)
}


#' Delete custom label(s)
#'
#' Removes the given custom label(s) from all genotypes
#' and deletes them from simulation memory.
#'
#' @param labels an vector containing the label numbers of the labels to be deleted
#'
#' @family label functions
#' @family deletor functions
#' @export
delete.label <- function(labels) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  labels <- convert.to.integer(labels,"Label identifiers")
  .Call(SXP_delete_label, sim.data$p, labels)
}


#' Delete marker effect set(s)
#'
#' Removes the given marker effect set(s) 
#' from simulation memory.
#'
#' @param effect.sets an vector containing the effect set identifiers of the 
#' effect sets to be deleted
#'
#' @family deletor functions
#' @export
delete.effects <- function(effect.sets) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  effect.sets <- convert.to.integer(effect.sets,"Effect set identifiers")
  .Call(SXP_delete_eff_set, sim.data$p, effect.sets)
}

#' OLD NAME | Delete marker effect set(s)
#'
#' ! This is the old name for \code{delete.effects}. It has been changed for consistency
#' with \code{load.effects}
#'
#' @seealso \link{delete.effects}
#' 
#' @keywords internal 
#' @export
delete.effect.set <- function(effect.sets) {
  return(delete.effects(effect.sets))
}


#' Delete recombination map(s)
#'
#' Removes the given recombination map(s)
#' from simulation memory.
#'
#' @param maps an vector containing the map identifiers of the maps to be deleted
#'
#' @family deletor functions
#' @export
delete.map <- function(maps) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  maps <- convert.to.integer(maps,"Map identifiers")
  .Call(SXP_delete_recombination_map, sim.data$p, maps)
}


#' OLD NAME | Delete recombination map(s)
#' 
#' ! This is the old name for \code{delete.map}. It has been changed for consistency
#' with \code{load.map}
#'
#' @seealso \link{delete.map}
#' 
#' @keywords internal 
#' @export
delete.recombination.map <- function(maps) {
  return(delete.map(maps))
}

