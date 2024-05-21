#' Clear the internal storage of all data
#'
#' \code{clear.simdata} frees and deletes all data stored in the
#' package's internal SimData struct (sim.data$p). Returns 0 on success.
#'
#' @family loader functions
#' @family deletor functions
#' @export
clear.simdata <- function() {
	if (is.null(sim.data$p)) { return(0L) }
	sim.data$p <- NULL
	gc()
	return(0L)
}


#' Delete a group's worth of genotypes
#'
#' \code{delete.group} finds the genotypes of the given group and removes
#' them from the SimData object's storage. A message is printed explaining
#' how many genotypes were deleted.
#' 
#' Multiple group input supported.
#'
#' @param groups an vector containing the group numbers of the groups to be deleted
#' @return 0 on success. An error is raised on failure.
#'
#' @family grouping functions
#' @family deletor functions
#' @export
delete.group <- function(groups) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_delete_group, sim.data$p, groups))
}


#' Delete custom label(s)
#'
#' \code{delete.label} removes the given custom label(s) from all genotypes
#' and deletes them from simulation memory.
#'
#' @param labels an vector containing the label numbers of the labels to be deleted
#' @return 0 on success. An error is raised on failure.
#'
#' @family label functions
#' @family deletor functions
#' @export
delete.label <- function(labels) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_delete_label, sim.data$p, labels))
}


#' Delete marker effect set(s)
#'
#' \code{delete.effect.set} removes the given marker effect set(s) 
#' from simulation memory.
#'
#' @param effect.sets an vector containing the effect set identifiers of the 
#' effect sets to be deleted
#' @return 0 on success. An error is raised on failure.
#'
#' @family deletor functions
#' @export
delete.effect.set <- function(effect.sets) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_delete_eff_set, sim.data$p, effect.sets))
}


#' Delete recombination map(s)
#'
#' \code{delete.recombination.map} removes the given recombination map(s)
#' from simulation memory.
#'
#' @param maps an vector containing the map identifiers of the maps to be deleted
#' @return 0 on success. An error is raised on failure.
#'
#' @family deletor functions
#' @export
delete.recombination.map <- function(maps) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_delete_recombination_map, sim.data$p, maps))
}
